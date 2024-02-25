// Copyright (c) 2023-2024 Steve Young
// Licensed under the MIT License

#include "SkeletonPointCollection.h"
#include "SuperPixel.h"
#include "SPixelData.h"
#include "Workspace.h"

SkeletonPointCollection::SkeletonPointCollection(SuperPixel* sp)
{
	skeleton_superpixel = sp;
	skeleton_point_set.clear();
}

SkeletonPointCollection::~SkeletonPointCollection()
{
	ClearPoints();
}

bool SkeletonPointCollection::FindSkeletonPoints()
{
	int i, j, k, l, rot, transition_count, corner_count;
	PointPair p;
	bool prev_state, local_bit, corner_set;
	unsigned char bit0, bit1, bit2;
	RectQuad boundingbox = skeleton_superpixel->GetWindow();
	int identifier = skeleton_superpixel->GetIdentifier();
	pixdata = skeleton_superpixel->GetPixelData();
	width = skeleton_superpixel->GetPixelData()->GetWidth();
	height = skeleton_superpixel->GetPixelData()->GetHeight();

	ClearPoints();
	for (j = boundingbox.y0; j <= boundingbox.y1; j++)
	{
		for (i = boundingbox.x0; i <= boundingbox.x1; i++)
		{
			if (identifier == pixdata->GetPixel(i, j))
			{
				transition_count = 0;
				corner_count = 0;
				corner_set = false;
				prev_state = true; // Start out with this as true, so first time checking rot=0 won't result in a false-to-true transition. 
				for (rot = 0; rot < 9; rot++)  // Rotate through top pixel clockwise back to the top pixel again.
				{
					bit0 = rot & 1;
					bit1 = (rot & 2) >> 1;
					bit2 = (rot & 4) >> 2;
					l = j + (bit2 * 2 - 1) * (-1) * ((bit1 & bit0) + bit1 - 1);  // sets l, relative to j, as -1 -1 0 1 1 1 0 -1.
					k = i + (bit2 * 2 - 1) * (-1) * (bit0 | bit1);  // sets k, relative to i, as 0 1 1 1 0 -1 -1 -1.
					local_bit = false;
					if ((k >= boundingbox.x0) && (k <= boundingbox.x1) && (l >= boundingbox.y0) && (l <= boundingbox.y1))
					{
						local_bit = (identifier == pixdata->GetPixel(k, l));
					}
					if (local_bit && (false == prev_state)) // This is a false-to-true transition.
					{
						transition_count++;
					}
					prev_state = local_bit;
					// Checking for case of a block of four pixels which would not otherwise get identified as intersection points.
					if (false == local_bit)  // Reset count to zero any time there is a not set bit.
					{
						corner_count = 0;
					}
					else {
						if (0 == bit0) // Top, right, bottom, or left.  Only start the count on one of these.
						{
							corner_count++;
						}
						else {  // A corner.
							if (corner_count > 0) // Only counts if the previous adjacent side was counted, so the counter is above zero.
							{
								corner_count++;
							}
						}
						if (corner_count > 2)
						{
							corner_set = true;
						}
					}
				}
				if (1 >= transition_count) // Endpoint.
				{
					p.x = i;
					p.y = j;
					AddSkeletonEndpoint(p);
				}
				else if ((2 == transition_count) && corner_set) // Rare instance of an intersetion point of four pixels forming a square.
				{
					p.x = i;
					p.y = j;
					AddSkeletonIntersection(p);
				}
				else if (transition_count > 2) // Intersection point.
				{
					p.x = i;
					p.y = j;
					AddSkeletonIntersection(p);
				}
			}
		}
	}
	if ((GetIntersections().size() == 0) && (GetEndpoints().size() == 0))
	{
		// If there are not intersections, add one.
		bool cont = true;
		for (j = boundingbox.y0; (j <= boundingbox.y1) && cont; j++)
		{
			for (i = boundingbox.x0; (i <= boundingbox.x1) && cont; i++)
			{
				if (identifier == pixdata->GetPixel(i, j))
				{
					p.x = i;
					p.y = j;
					AddSkeletonIntersection(p);
					cont = false;
				}
			}
		}

	}
	return true;
}

bool SkeletonPointCollection::FindSkeletonLinks()
{

	int identifier = skeleton_superpixel->GetIdentifier();
	int neighbors, core_neighbors, exclude;
	int pos, new_pos, core_pos;
	Skeleton_Link* link = NULL;
	PointPair pixel;
	PointPair core_pixel;
	std::set<SkeletonPointIndex>::iterator it;
	std::set<int>::iterator adjacent_iterator;

	int width = pixdata->GetWidth();
	Skeleton_Point* current_point;
	Skeleton_Point* core_point;
	std::set<int> seen_points;
	std::set<int> intersections = GetIntersections();

	//for (it = endpoints.begin(); it != endpoints.end(); ++it)
	//{
	//	bool seen = false;  // Has this point already shown up in an earlier link?
	//	pos = it->point; // The integer position of the point.
	//	current_point = GetByIdentifier(pos);
	//	if (NULL == current_point)
	//	{
	//		throw std::runtime_error("Skeleton endpoint not included in index.\n");
	//		return false;
	//	}
	//	for (int i = 0; i < 8; i++)
	//	{
	//		if (NULL != current_point->dir[i])
	//		{
	//			seen = true;
	//		}
	//	}
	//	if (false == seen) // For endpoints, no need to do further examination if it has already been made part of a link.
	//	{
	//		link = new Skeleton_Link();
	//		if (NULL == link)
	//		{
	//			throw std::runtime_error("Failed to create new skeleton link object.\n");
	//			return false;
	//		}
	//		link->end1 = pos; // Initially, assign iterator ('it') to end1, but may need to reverse these later, so end1 < end2 in all cases (to avoid duplication).
	//		link->link_points.insert(pos);
	//		link->line.push_back(pos);
	//		link->distance = 0;
	//		bool found = true;
	//		bool first = true;
	//		while (found)
	//		{
	//			found = false;
	//			pixel = skeleton_superpixel->Pos2XY(pos);
	//			neighbors = pixeldata->GetNeighbors(identifier, pixel);
	//			// Neighbor values around point o:
	//			//
	//			// 128   1   2
	//			//  64   o   4
	//			//  32  16   8
	//			//
	//			// First, check to see whether any adjacent points are part of the set.  If so, that will be the next point.
	//			for (int i = 0; i < 4; i++)
	//			{
	//				int j = 1 << (i * 2); // Check 1, 4, 16, 64
	//				if (j == (neighbors & j))
	//				{
	//					new_pos = pos + (i & 1) * (2 - i) + (1 - (i & 1)) * (i - 1) * width; // Should convert 0-4 into the new pos for the neighbor.
	//					if (link->link_points.end() == link->link_points.find(new_pos)) // Test for whether this new point has already been encountered in this link.
	//					{
	//						link->link_points.insert(new_pos);
	//						link->line.push_back(new_pos);
	//						if (first)
	//						{
	//							current_point->dir[i * 2] = link; // Set the directional link for the existing point.
	//						}
	//						link->distance += 1.0; // Distance increases by 1.0 for adjancent pixels.
	//						found = true;
	//						break;
	//					}
	//				}
	//			}
	//			if (false == found) // None of the neighbors were the adjacent pixels.  Now look at diagonals.
	//			{
	//				for (int i = 0; i < 4; i++)
	//				{
	//					int j = 1 << (1 + i * 2); // Check 2, 8, 32, 128
	//					if (j == (neighbors & j))
	//					{
	//						new_pos = pos + (1 - (i & 2)) + (2 * ((i & 1) ^ ((i & 2) >> 1)) - 1) * width; // Should convert 0-4 into new pos.
	//						if (link->link_points.end() == link->link_points.find(new_pos)) // Test for whether this new point has already been encountered in this link.
	//						{
	//							link->link_points.insert(new_pos);
	//							link->line.push_back(new_pos);
	//							if (first)
	//							{
	//								current_point->dir[i * 2 + 1] = link;
	//							}
	//							link->distance += 1.4142;  // Distance increases by the hypotenuse for diagonal points.
	//							found = true;
	//							break;
	//						}
	//					}
	//				}
	//			}
	//			if (found)  // New pixel has been found as part of the link.
	//			{
	//				pos = new_pos;
	//				// Need to test to see if this is an endpoint or intersection.  If so, no more tracing for this link.
	//				if (TestEndpointOrIntersection(pos))
	//				{
	//					found = false; // To end the while loop.
	//					link->end2 = pos;
	//					// Put directional data into the endpoint or intersection.
	//					current_point = GetByIdentifier(pos);
	//					// Set the direction pointers for the current_point.
	//					int line_size = link->line.size();
	//					if (line_size > 1)
	//					{
	//						int j = GetDirection(link->line[line_size - 1], link->line[line_size - 2], width);
	//						if (j >= 0)
	//						{
	//							current_point->dir[j] = link;
	//						}
	//						else {
	//							throw std::runtime_error("Non-contiguous final points for skeleton line.\n");
	//							return false;
	//						}
	//						if (line_size > 2) // Sometimes the last two points on a line will touch the endpoint or intersection.
	//						{
	//							j = GetDirection(link->line[line_size - 1], link->line[line_size - 3], width);
	//							if (j >= 0)
	//							{
	//								current_point->dir[j] = link;
	//							}
	//						}
	//					}
	//					else {
	//						throw std::runtime_error("Too few points in skeleton link.\n");
	//						return false;
	//					}
	//					if (link->end1 > link->end2) // Need to reverse these for canonical order.
	//					{
	//						link->end2 = link->end1;
	//						link->end1 = pos;
	//						// Need to reverse order of line vector.
	//						std::reverse(link->line.begin(), link->line.end());
	//					}
	//				}
	//			}
	//			else if (first) // Single pixel makes up this skeleton.
	//			{
	//				link->end2 = pos;
	//				links.insert(*link);
	//			}
	//			else {
	//				throw std::runtime_error("Skeleton line ends without endpoint or intersection.\n");
	//				return false;
	//			}
	//			first = false;
	//		}
	//	}
	//}
	for (it = skeleton_point_set.begin(); it != skeleton_point_set.end(); ++it)
	{
		seen_points.clear();
		core_pos = it->point; // The integer position of the point.
		core_point = it->sp;
		if (NULL == core_point)
		{
			throw std::runtime_error("Skeleton point not included in index.\n");
			return false;
		}
		std::set<int> adjacent_intersections = GetAdjacentIntersections(core_pos);
		int num_adjacent_intersections = adjacent_intersections.size();

		for (int i = 0; i < 8; i++)
		{
			if (NULL != core_point->dir[i])
			{
				seen_points.insert(core_point->dir[i]->link_points.begin(), core_point->dir[i]->link_points.end());
			}
		}

		core_pixel = skeleton_superpixel->Pos2XY(core_pos);
		core_neighbors = pixdata->GetNeighbors(identifier, core_pixel);
		exclude = 0; // No exclusions to start with.  This is for use to avoid a turn-back after the first point of a link.
		for (int ii = 0; ii < 8; ii++)
		{
			std::set<int>::iterator intersection_iterator;
			for (intersection_iterator = intersections.begin(); intersection_iterator != intersections.end(); ++intersection_iterator)
			{
				int int_point = *intersection_iterator;
				seen_points.erase(int_point);  // Remove all intersection points, to allow for multiple links from this core_point to the same intersection.
			}
			int pix_i = ii * 2 - 7 * ((ii & 4) >> 2);  // Takes pix_i through the adjacent pixels (0, 2, 4, 6) then diagonals (1, 3, 5, 7).
			if (NULL == core_point->dir[pix_i])  // We don't want to repeat a link that has already been done.
			{
				int jj = 1 << pix_i;  // Neighbors bitmask.
				// pix_i bit2 bit1 bit0 dx dy
				// 0     0    0    0    0  -w
				// 1     0    0    1    1  -w
				// 2     0    1    0    1   0
				// 3     0    1    1    1   w
				// 4     1    0    0    0   w
				// 5     1    0    1    -1  w
				// 6     1    1    0    -1  0
				// 7     1    1    1    -1  -w

				// dx = (bit1|bit0)*(1-2*bit2)
				// dy = w*(bit0|(~bit1))*(-1 + 2*(bit2^bit1))
				unsigned char bit0, bit1, bit2;
				bit0 = pix_i & 1;
				bit1 = (pix_i & 2) >> 1;
				bit2 = (pix_i & 4) >> 2;
				bool skip_adjacent = false;  // Flag for whether this pixel shoud be skipped because it is a diagonal adjacent to an adjacent intersection.
				if ((1 == bit0) && (num_adjacent_intersections > 0)) // Diagonal and there are adjacent intersections.
				{
					new_pos = core_pos + (bit1 | bit0) * (1 - 2 * bit2) + (bit0 | (~bit1) & 1) * (-1 + 2 * ((bit2 ^ bit1) & 1)) * width; // Should convert pix_i to the new pos for the neighbor.

					for (adjacent_iterator = adjacent_intersections.begin(); adjacent_iterator != adjacent_intersections.end(); ++adjacent_iterator)
					{
						int diff = abs(new_pos - *adjacent_iterator);
						if ((width == diff) || (1 == diff)) // Will be true if new_pos is adjacent to the point that is adjacent to the core_pos.
						{
							skip_adjacent = true;
							break;
						}
					}
				}
				if (false == skip_adjacent) { // Make sure diagonal is not adjacent to an adjacent intersection.
					if (jj == (core_neighbors & jj)) // There is a connecting pixel.
					{
						link = new Skeleton_Link();
						if (NULL == link)
						{
							throw std::runtime_error("Failed to create new skeleton link object.\n");
							return false;
						}
						link->end1 = core_pos; // Initially, assign core_pos to end1, but may need to reverse these later, so end1 < end2 in all cases (to avoid duplication).
						link->link_points.insert(core_pos);  // This is the 'core' pixel, not the first connecting one.
						link->line.push_back(core_pos);
						seen_points.insert(core_pos);
						link->distance = 0;

						bool found = true;
						bool first = true;
						pos = core_pos;
						while (found)
						{
							found = false;
							pixel = skeleton_superpixel->Pos2XY(pos);
							if (first) // For the first round, use the pixel identified above using the pix_i value.
							{
								found = true;

								new_pos = pos + (bit1 | bit0) * (1 - 2 * bit2) + (bit0 | (~bit1) & 1) * (-1 + 2 * ((bit2 ^ bit1) & 1)) * width; // Should convert pix_i to the new pos for the neighbor.

								link->link_points.insert(new_pos);
								link->line.push_back(new_pos);
								seen_points.insert(new_pos);
								core_point->dir[pix_i] = link; // Set the directional link for the first point.
								link->distance += 1.0 + 0.4142 * (pix_i & 1); // Distance increases by 1.0 for adjancent pixels, plus more for diagonal.
								if (0 == (bit0 & 1)) // Not a diagonal.
								{
									exclude = (40 + bit1 * 120 + bit2 * 90) % 240;
								}
								else {
									exclude = 0;
								}
							}
							else {
								neighbors = pixdata->GetNeighbors(identifier, pixel);
								// First, check to see whether any adjacent points are part of the set.  If so, that will be the next point.
								for (int i = 0; i < 4; i++)
								{
									int j = 1 << (i * 2); // Check 1, 4, 16, 64
									if (j == (neighbors & j))
									{
										new_pos = pos + (i & 1) * (2 - i) + (1 - (i & 1)) * (i - 1) * width; // Should convert 0-4 into the new pos for the neighbor.
										if (seen_points.end() == seen_points.find(new_pos)) // Test whether this point has been seen for this intersection.
										{
											link->link_points.insert(new_pos);
											link->line.push_back(new_pos);
											seen_points.insert(new_pos);
											link->distance += 1.0; // Distance increases by 1.0 for adjancent pixels.
											found = true;
											if (3 == link->line.size())
											{
												int dir_index = GetDirection(core_pos, link->line[link->line.size() - 1], width);
												if (dir_index >= 0)
												{
													if (NULL == core_point->dir[dir_index])
													{
														core_point->dir[dir_index] = link;
													}
												}
											}
											if (4 == link->line.size())
											{
												seen_points.erase(core_pos); // To allow for a loop around.
											}
											break;
										}
									}
								}
								if (false == found) // None of the neighbors were the adjacent pixels.  Now look at diagonals.
								{
									for (int i = 0; i < 4; i++)
									{
										int j = 1 << (1 + i * 2); // Check 2, 8, 32, 128
										if (j == (neighbors & j) && ((exclude == 0) || (j != (exclude & j))))
										{
											new_pos = pos + (1 - (i & 2)) + (2 * ((i & 1) ^ ((i & 2) >> 1)) - 1) * width; // Should convert 0-4 into new pos.

											if (seen_points.end() == seen_points.find(new_pos)) // Test whether this point has been seen for this intersection.
											{
												link->link_points.insert(new_pos);
												link->line.push_back(new_pos);
												seen_points.insert(new_pos);
												link->distance += 1.4142;  // Distance increases by the hypotenuse for diagonal points.
												found = true;
												if (3 == link->line.size())
												{
													int dir_index = GetDirection(core_pos, link->line[link->line.size() - 1], width);
													if (dir_index >= 0)
													{
														if (NULL == core_point->dir[dir_index])
														{
															core_point->dir[dir_index] = link;
														}
													}
												}
												if (4 == link->line.size())
												{
													seen_points.erase(core_pos); // To allow for a loop around.
												}
												break;
											}
										}
									}
								}
								exclude = 0;
							}
							if (found)  // New pixel has been found as part of the link.
							{
								pos = new_pos;
								// Need to test to see if this is an endpoint or intersection.  If so, no more tracing for this link.
								if (TestEndpointOrIntersection(pos))
								{
									found = false; // To end the while loop.
									link->end2 = pos;
									// Put directional data into the endpoint or intersection.
									current_point = GetSPByIdentifier(pos);

									// Set the direction pointers for the current_point.
									int line_size = link->line.size();
									if (line_size > 1)
									{
										int j = GetDirection(link->line[line_size - 1], link->line[line_size - 2], width);
										if (j >= 0)
										{
											current_point->dir[j] = link;
										}
										else {
											throw std::runtime_error("Non-contiguous final points for skeleton line.\n");
											return false;
										}
										if (line_size > 2) // Sometimes the last two points on a line will touch the endpoint or intersection.
										{
											j = GetDirection(link->line[line_size - 1], link->line[line_size - 3], width);
											if (j >= 0)
											{
												current_point->dir[j] = link;
											}
										}
									}
									else {
										throw std::runtime_error("Too few points in skeleton link.\n");
										return false;
									}
									if (link->end1 > link->end2) // Need to reverse these for canonical order.
									{
										link->end2 = link->end1;
										link->end1 = pos;
										// Need to reverse order of line vector.
										std::reverse(link->line.begin(), link->line.end());
									}
								}
							}
							else {
								throw std::runtime_error("Skeleton line from intersection ends without endpoint or intersection.\n");
								return false;
							}
							first = false;
						}
					}
				}
			}
		}
	}
	return true;
}

bool SkeletonPointCollection::FindLongestPath()
{
	std::set<SkeletonPointIndex>::iterator point_it;
	std::set<SkeletonPointIndex>::iterator far_point_it;
	std::set<int>::iterator leads_it;
	Skeleton_Point* other_end;
	Skeleton_Point* local_lead;
	Skeleton_Link* local_link;
	Skeleton_Point* sp;
	Skeleton_Point* far_sp;
	SPixelData* pixeldata = skeleton_superpixel->GetWorkspace()->GetPixeldata();
	int local_point;
	std::set<int> leads;  // A set of leads going in different directions.
	std::vector<int> new_leads;  // New leads to be followed after processing existing ones.
	std::vector<int> long_path;
	std::set<Skeleton_Link*> link_path; // Temporary list of links making up the long_path;
	float long_path_dist;
	float longest;
	bool done = false;
	bool first = true;

	ClearLongPathFlags();
	longest = 0.0;
	longest_path.clear();
	first = true;

	if (1 == skeleton_point_set.size())
	{
		bool has_link = false;
		for (int i = 0; (i < 8)&&(!has_link); ++i)
		{
			point_it = skeleton_point_set.begin();
			if (NULL != point_it->sp->dir[i])
			{
				has_link = true;
			}
		}
		if (!has_link) // Only one point, no links to it.
		{
			longest_path.clear();
			local_point = skeleton_point_set.begin()->point;
			longest_path.push_back(local_point);
			longest_paths.clear();
			longest_paths.insert(longest_path);
			return true;
		}
	}

	while (false == done)
	{
		long_path_dist = 0.0;
		long_path.clear();

		for (point_it = skeleton_point_set.begin(); point_it != skeleton_point_set.end(); ++point_it)
		{
			// First part is to popluate all temp_dist values with shortest distances from point_it->sp.
			ResetTempDists();  // Set all temp_dist values to -1 (so we know they have not been visted yet).
			sp = point_it->sp;
			sp->temp_dist = 0.0;
			leads.clear();
			new_leads.clear();
			new_leads.push_back(sp->point);

			while (new_leads.size() > 0)  // Keep looping through until no more leads are being added.
			{
				leads.clear();
				leads.insert(new_leads.begin(), new_leads.end());
				new_leads.clear();
				for (leads_it = leads.begin(); leads_it != leads.end(); ++leads_it)
				{
					local_point = *leads_it;
					local_lead = GetSPByIdentifier(local_point);
					for (int i = 0; i < 8; ++i)
					{
						local_link = local_lead->dir[i];
						if ((NULL != local_link) && (false == local_link->used_in_long_path))
						{
							if (local_link->end1 == local_point)
							{
								other_end = GetSPByIdentifier(local_link->end2);
							}
							else {
								other_end = GetSPByIdentifier(local_link->end1);
							}
							if (other_end->temp_dist < 0)
							{
								other_end->temp_dist = local_lead->temp_dist + local_link->distance;
								new_leads.push_back(other_end->point);
							}
							else if ((local_lead->temp_dist + local_link->distance) < other_end->temp_dist)
							{
								other_end->temp_dist = local_lead->temp_dist + local_link->distance;
								new_leads.push_back(other_end->point);
							}
						}
					}
				}
			}
			// Now go through all points to find the one with the longest distance.
			float local_distance = 0.0;
			int local_end = sp->point;
			for (far_point_it = skeleton_point_set.begin(); far_point_it != skeleton_point_set.end(); ++far_point_it)
			{
				far_sp = far_point_it->sp;
				if (far_sp->temp_dist > local_distance)
				{
					local_distance = far_sp->temp_dist;
					local_end = far_sp->point;
				}
			}
			far_sp = GetSPByIdentifier(local_end);

			if (local_distance > long_path_dist)  // No need to work out the exact path unless this is longer than others seen so far.
			{
				// far_sp will move from the ultimate terminal point to the beginning point.
				// local_end is the position of far_sp.
				// sp is the beginning point of the path (which is arrived at last, since we are backtracking from far_sp).
				// local_lead is the next point along the path backtracking towards sp.
				// local_lead_pos is the position of local_lead.

				int local_lead_pos;
				long_path.clear();
				link_path.clear();  // Reset link_path.
				long_path.push_back(local_end);
				long_path_dist = local_distance;
				while (local_end != sp->point)
				{
					int i = FindDirSmallestTempDist(far_sp);
					if (i < 0)
					{
						throw std::runtime_error("Error backtracking longest path for skeleton.\n");
						return false;
					}
					bool reverse = false;  // True if tracing back requires going through line from reverse order.
					link_path.insert(far_sp->dir[i]);  // Add link to the link_path.
					if (far_sp->dir[i]->end1 == local_end)  // Need to look to end2.
					{
						local_lead_pos = far_sp->dir[i]->end2;
					}
					else { // Need to look to end1.
						local_lead_pos = far_sp->dir[i]->end1;
						reverse = true;
					}
					local_lead = GetSPByIdentifier(local_lead_pos);
					int line_index = 1;  // Not zero because the point at location zero (or size-1 if reverse is true) is already in long_path.
					int line_size = far_sp->dir[i]->line.size();
					if (line_size > 2)
					{
						if (reverse)
						{
							line_index = line_size - 2;
						}
						while (long_path[0] != local_lead_pos)
						{
							long_path.insert(long_path.begin(), far_sp->dir[i]->line[line_index]);
							if (reverse)
							{
								line_index--;
							}
							else {
								line_index++;
							}
							if ((line_index > line_size) || (line_index < -1))
							{
								throw std::runtime_error("Failed to find next point in skeleton line segment.\n");
								return false;
							}
						}
					}
					else {
						long_path.insert(long_path.begin(), local_lead_pos);
					}
					local_distance = local_lead->temp_dist;
					local_end = local_lead_pos;
					far_sp = GetSPByIdentifier(local_end);
				}
			}
		}
		if (first || (long_path_dist > 5.0))
		{
			if (first)
			{
				longest = long_path_dist;
				longest_path = long_path;
				first = false;
			}
			longest_paths.insert(long_path);
			int prev = 0;
			for (std::vector<int>::iterator it = long_path.begin(); it != long_path.end(); ++it)
			{
				if (long_path.begin() == it)
				{
					prev = *it;
				}
				else {
					sp = GetSPByIdentifier(*it);
					if (NULL != sp)
					{
						for (int i = 0; i < 8; ++i)
						{
							if ((NULL != sp->dir[i]) && (false == sp->dir[i]->used_in_long_path))
							{
								if (((sp->dir[i]->end1 == prev) && (sp->dir[i]->end2 == *it)) ||
									((sp->dir[i]->end2 == prev) && (sp->dir[i]->end1 == *it)))
								{
									std::set<Skeleton_Link*>::iterator find_it;
									find_it = link_path.find(sp->dir[i]);
									if (find_it != link_path.end())
									{
										sp->dir[i]->used_in_long_path = true;
									}
								}
							}
						}
						prev = *it;
					}
				}
			}
		}
		else {
			done = true;
			// Still need to pick up orphaned links (which may be pointing to the same point on both ends.
			for (point_it = skeleton_point_set.begin(); point_it != skeleton_point_set.end(); ++point_it)
			{
				sp = point_it->sp;
				for (int i = 0; i < 8; ++i)
				{
					local_link = sp->dir[i];
					if ((NULL != local_link) && (false == local_link->used_in_long_path) && (local_link->distance > 5.0))
					{
						long_path.clear();
						for (std::vector<int>::iterator it = local_link->line.begin(); it != local_link->line.end(); ++it)
						{
							long_path.push_back(*it);
						}
						longest_paths.insert(long_path);
						local_link->used_in_long_path = true;
					}
				}
			}
		}
	}
	return true;
}

bool SkeletonPointCollection::AddSkeletonEndpoint(PointPair p)
{
	Skeleton_Point* new_point = new Skeleton_Point();
	if (NULL == new_point)
	{
		throw std::runtime_error("Failed to allocate memory for new skeleton endpoint.\n");
		return false;
	}
	new_point->point = skeleton_superpixel->XY2Pos(p);
	for (int i = 0; i < 8; i++)
	{
		new_point->dir[i] = NULL;
	}
	new_point->endpoint = true;
	InsertPointIntoIndex(new_point);
	return true;
}

bool SkeletonPointCollection::AddSkeletonIntersection(PointPair p)
{
	Skeleton_Point* new_point = new Skeleton_Point();
	if (NULL == new_point)
	{
		throw std::runtime_error("Failed to allocate memory for new skeleton intersection point.\n");
		return false;
	}
	new_point->point = skeleton_superpixel->XY2Pos(p);
	for (int i = 0; i < 8; i++)
	{
		new_point->dir[i] = NULL;
	}
	new_point->endpoint = false;
	InsertPointIntoIndex(new_point);
	return true;
}

bool SkeletonPointCollection::RecalcDistances()
{
	std::set<SkeletonPointIndex>::iterator it;
	for (it = skeleton_point_set.begin(); it != skeleton_point_set.end(); ++it)
	{
		Skeleton_Point* local_point = it->sp;
		int local_pos = local_point->point;
		for (int i = 0; i < 8; ++i)
		{
			if (NULL != local_point->dir[i])
			{
				Skeleton_Link* link = local_point->dir[i];
				if (local_pos == link->end1) // Only recalculate each link once.
				{
					float dist = 0; // Distance for current point.
					float dist_1 = -1; // Distance for current point minus 1.
					float dist_2 = -1; // Distance for current point minus 2.
					int size = link->line.size();
					int diff;
					for (int index = 0; index < size; ++index)
					{
						if (dist_1 >= 0)
						{
							diff = abs(link->line[index] - link->line[index - 1]);
							if ((diff == 1) || (diff == width)) // Previous point is adjacent to current point.
							{
								dist = dist + 1.0;;
							}
							else { // Previous point is diagonal to current point.
								dist = dist + 1.4142;
							}
							if (dist_2 >= 0)
							{
								diff = abs(link->line[index] - link->line[index - 2]);
								if ((diff == 1) || (diff == width)) // Two points back was adjacent to current point.  Can this happen?
								{
									if (dist > (dist_2 + 1.0))
									{
										dist = dist_2 + 1.0;
									}
								}
								else if ((diff == (width + 1)) || (diff == (width - 1))) // Two points back was diagonal to current point.
								{
									if (dist > (dist_2 + 1.4142))
									{
										dist = dist_2 + 1.4142;
									}
								}
							}
						}
						dist_2 = dist_1; // Prepare for next iteration.
						dist_1 = dist;
					}
					link->distance = dist;
				}
			}
		}
	}
	return true;
}

int SkeletonPointCollection::FindDirSmallestTempDist(Skeleton_Point* sp)
{
	// Using sp, look at all non-NULL links to it, which have not already been used in a long path, and find the one that leads to the smallest temp_dist.
	int ret = -1;
	float smallest = -1.0;
	int pos = sp->point;
	int other_pos = -1;
	int ret_pos = -1;
	Skeleton_Point* other_sp = NULL;
	for (int i = 0; i < 8; ++i)
	{
		if ((NULL != sp->dir[i]) && (false == sp->dir[i]->used_in_long_path))
		{
			if (sp->dir[i]->end1 == pos)
			{
				other_pos = sp->dir[i]->end2;
			}
			else {
				other_pos = sp->dir[i]->end1;
			}
			other_sp = GetSPByIdentifier(other_pos);
			if ((smallest < 0.0) || (other_sp->temp_dist < smallest))
			{
				smallest = other_sp->temp_dist;
				ret = i;
				ret_pos = other_pos;
			}
		}
	}
	if (ret_pos > -1)  // We know that ret_pos is the point the link will go to, but we may not yet know the shortest link to it.
	{
		smallest = -1.0;
		for (int i = 0; i < 8; ++i)
		{
			if ((NULL != sp->dir[i]) && (false == sp->dir[i]->used_in_long_path))
			{
				if (sp->dir[i]->end1 == pos)
				{
					other_pos = sp->dir[i]->end2;
				}
				else {
					other_pos = sp->dir[i]->end1;
				}
				if ((other_pos == ret_pos) && ((smallest < 0) || (smallest > sp->dir[i]->distance)))
				{
					smallest = sp->dir[i]->distance;
					ret = i;
				}
			}
		}
	}

	return ret;
}

std::set<int> SkeletonPointCollection::GetEndpoints()
{
	std::set<int> endpoint_set;
	endpoint_set.clear();
	std::set<SkeletonPointIndex>::iterator it;
	for (it = skeleton_point_set.begin(); it != skeleton_point_set.end(); ++it)
	{
		if (it->sp->endpoint)
		{
			endpoint_set.insert(it->point);
		}
	}
	return endpoint_set;
}

std::set<int> SkeletonPointCollection::GetIntersections()
{
	std::set<int> intersection_set;
	intersection_set.clear();
	std::set<SkeletonPointIndex>::iterator it;
	for (it = skeleton_point_set.begin(); it != skeleton_point_set.end(); ++it)
	{
		if (false == it->sp->endpoint)
		{
			intersection_set.insert(it->point);
		}
	}
	return intersection_set;
}

std::set<int> SkeletonPointCollection::GetAdjacentIntersections(int pos)
{
	std::set<int> ret;
	ret.clear();

	int y = (int)(pos / width);
	int x = pos - (y * width);
	int identifier = pixdata->GetPixel(x, y);

	if (y > 0)
	{
		if (identifier == pixdata->GetPixel(x, y - 1))
		{
			ret.insert(pos - width);
		}
	}
	if (y < (height - 1))
	{
		if (identifier == pixdata->GetPixel(x, y + 1))
		{
			ret.insert(pos + width);
		}
	}
	if (x > 0)
	{
		if (identifier == pixdata->GetPixel(x - 1, y))
		{
			ret.insert(pos - 1);
		}
	}
	if (x < (width - 1))
	{
		if (identifier == pixdata->GetPixel(x + 1, y))
		{
			ret.insert(pos + 1);
		}
	}
	return ret;
}

bool SkeletonPointCollection::TestEndpointOrIntersection(int pos)
{
	Skeleton_Point* point = GetSPByIdentifier(pos);
	if (NULL != point)
	{
		return true;
	}
	return false;
}

Skeleton_Point* SkeletonPointCollection::GetSPByIdentifier(int key)
{

	Skeleton_Point* ret = NULL;
	SkeletonPointIndex index;
	index.point = key;
	std::set<SkeletonPointIndex>::iterator it;
	it = skeleton_point_set.find(index);
	if (it != skeleton_point_set.end())
	{
		index = *it;
		ret = index.sp;
	}
	return ret;
}

bool SkeletonPointCollection::InsertPointIntoIndex(Skeleton_Point* sp)
{
	SkeletonPointIndex index;
	index.point = sp->point;
	index.sp = sp;
	skeleton_point_set.insert(index);
	return true;
}

int SkeletonPointCollection::GetDirection(int p1, int p2, int width)
{
	// p1 is the main point.
	// p2 is a point contiguous with p1.  The direction of p2 is 0 if up, and moves through 7 clockwise around p1.
	// if not contiguous, return is -1.

	int ret = -1;
	int diff = p1 - p2;
	if (width == diff)
	{
		ret = 0;
	}
	else if ((width - 1) == diff)
	{
		ret = 1;
	}
	else if (-1 == diff)
	{
		ret = 2;
	}
	else if ((-width - 1) == diff)
	{
		ret = 3;
	}
	else if (-width == diff)
	{
		ret = 4;
	}
	else if ((1 - width) == diff)
	{
		ret = 5;
	}
	else if (1 == diff)
	{
		ret = 6;
	}
	else if ((width + 1) == diff)
	{
		ret = 7;
	}

	return ret;
}

std::vector<int>* SkeletonPointCollection::GetLongestPath()
{
	return &longest_path;
}

std::set<std::vector<int>> SkeletonPointCollection::GetLongestPaths()
{
	return longest_paths;
}

bool SkeletonPointCollection::ClearPoints()
{
	skeleton_point_set.clear();
	return true;
}

bool SkeletonPointCollection::ResetTempDists()
{
	std::set<SkeletonPointIndex>::iterator it;
	for (it = skeleton_point_set.begin(); it != skeleton_point_set.end(); ++it)
	{
		SkeletonPointIndex spi = *it;
		spi.sp->temp_dist = -1.0;
	}
	return true;
}

bool SkeletonPointCollection::ClearLongPathFlags()
{
	std::set<SkeletonPointIndex>::iterator it;
	for (it = skeleton_point_set.begin(); it != skeleton_point_set.end(); ++it)
	{
		for (int i = 0; i < 8; ++i)
		{
			if (NULL != it->sp->dir[i])
			{
				it->sp->dir[i]->used_in_long_path = false;
			}
		}
	}
	return true;
}
