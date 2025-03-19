// Copyright (c) 2023-2025 Steve Young
// Licensed under the MIT License

#include "SuperPixel.h"
#include "Cell.h"
#include "Workspace.h"
#include "SkeletonPointCollection.h"

SuperPixel::SuperPixel(int id, GradData* graddat, SPixelData* pixdat, PointPair point, SuperPixel* n, SuperPixel* p, WorkSpace* ws, SuperPixelType t)
{
	seed = point;
	boundingbox.x0 = seed.x;
	boundingbox.y0 = seed.y;
	boundingbox.x1 = seed.x;
	boundingbox.y1 = seed.y;
	size = 1;
	next = n;
	prev = p;
	level_complete = 0;
	prevsize = 1;
	identifier = id;
	pixeldata = pixdat;
	gradientdata = graddat;
	workspace = ws;
	image_width = graddat->GetWidth();
	image_height = graddat->GetHeight();
	if (NULL != prev)
	{
		prev->SetNext(this);
	}
	if (NULL != next)
	{
		next->SetPrev(this);
	}
	pixeldata->SetPixel(seed.x, seed.y, id);
	EdgePixels.clear();
	EdgePixels.insert(XY2Pos(seed));
	Neighbors.clear();
	Vertices.clear();
	AveError = 0;
	num_paths = 0;
	type = t;
	EdgePixelsCurrent = false;
	fill_image.clear();
}

SuperPixel::SuperPixel(const SuperPixel& tsp, GradData* graddat, SPixelData* pixdat, SuperPixel* n, SuperPixel* p, WorkSpace* ws, SuperPixelType t)
{
	seed = tsp.seed;
	boundingbox = tsp.boundingbox;
	size = tsp.size;
	next = n;
	prev = p;
	level_complete = 0;
	prevsize = 1;
	identifier = tsp.identifier;
	pixeldata = pixdat;
	gradientdata = graddat;
	workspace = ws;
	image_width = tsp.image_width;
	image_height = tsp.image_height;
	if (NULL != prev)
	{
		prev->SetNext(this);
	}
	if (NULL != next)
	{
		next->SetPrev(this);
	}
	EdgePixels.clear();
	EdgePixels.insert(XY2Pos(seed));
	Neighbors.clear();
	Vertices.clear();
	AveError = 0;
	AveColor = tsp.AveColor;
	num_paths = 0;
	type = t;
	EdgePixelsCurrent = false;
	fill_image.clear();
}

SuperPixel::SuperPixel(int id, SPixelData* pixdat, PointPair point, SuperPixel* n, SuperPixel* p, WorkSpace* ws, SuperPixelType t)
{
	seed = point;
	boundingbox.x0 = seed.x;
	boundingbox.y0 = seed.y;
	boundingbox.x1 = seed.x;
	boundingbox.y1 = seed.y;
	size = 1;
	next = n;
	prev = p;
	level_complete = 0;
	prevsize = 1;
	identifier = id;
	pixeldata = pixdat;
	gradientdata = NULL;
	workspace = ws;
	image_width = pixdat->GetWidth();
	image_height = pixdat->GetHeight();
	if (NULL != prev)
	{
		prev->SetNext(this);
	}
	if (NULL != next)
	{
		next->SetPrev(this);
	}
	pixeldata->SetPixel(seed.x, seed.y, id);
	EdgePixels.clear();
	EdgePixels.insert(XY2Pos(seed));
	Neighbors.clear();
	Vertices.clear();
	AveError = 0;
	num_paths = 0;
	type = t;
	EdgePixelsCurrent = false;
	fill_image.clear();
}

SuperPixel::~SuperPixel()
{
	if (NULL != prev)
	{
		prev->SetNext(next);
	}
	if (NULL != next)
	{
		next->SetPrev(prev);
	}
	EdgePixels.clear();
	Vertices.clear();
	for (std::set<int>::iterator iter = Neighbors.begin(); iter != Neighbors.end(); ++iter)
	{
		SuperPixel* nsp = GetByIdentifier(*iter);
		if (NULL != nsp)
		{
			nsp->SetNeighbors();  // Assumes that instances of "identifier" in the pixeldata have been removed already.
		}
	}
	Neighbors.clear();
	if (NULL != skeleton_collection)
	{
		delete skeleton_collection;
	}
}

int SuperPixel::Grow(unsigned char value, bool limit, bool mode, RectQuad box, SPixelData* mask, int mask_value, bool diagonals)
// Return identifier if there is more growing to be done at this value.
// Return 0 if there is no more growing to be done at this value.
// Return the identifier associated with another SuperPixel if it is encountered.
// Return -1 on error.
// Mode is: ordinary growing = true, return when encountering another SuperPixel = false
// Limit is: limited to bounding box (if mask is NULL).  If mask is not NULL, use it (with mask_value).
// Diagonals is: can the set grow to diagonal pixels in addition to vertically and horizontally.
//
// Modified from original in that it does not grow diagonally.
{
	if (SPType_Plain != type)
	{
		throw std::runtime_error("Attempting to grow non-plain SuperPixel.\n");
		return -1;
	}
	bool use_mask = false;
	int ret = 0;
	if ((level_complete > 0) && (value <= level_complete))
	{
		return 0;
	}
	if (mask != NULL)
	{
		use_mask = true;
		limit = false;
	}
	if (false == limit)
	{
		box.x0 = 0;
		box.y0 = 0;
		box.x1 = pixeldata->GetWidth() - 1;
		box.y1 = pixeldata->GetHeight() - 1;
	}
	if (false == EdgePixelsCurrent)
	{
		FindEdgePixels();
	}
	if (EdgePixels.size() > size)
	{
		throw std::runtime_error("Edgepixel number exceeds size of SuperPixel.\n");
	}
	if (EdgePixels.size() > 0)
	{
		std::vector<int> NewEdgePixels;  // Measured 36% faster on test image using vector instead of set.
		NewEdgePixels.clear();
		std::set<int>::iterator it;
		for (it = EdgePixels.begin(); it != EdgePixels.end(); ++it)
		{
			int pos = *it;
			PointPair xy = Pos2XY(pos);
			int i = xy.x;
			int j = xy.y;
			bool surrounded = true;
			for (int l = j - 1; l <= j + 1; l++)
			{
				for (int k = i - 1; k <= i + 1; k++)
				{
					if ((l >= box.y0) && (k >= box.x0) && (l <= box.y1) && (k <= box.x1) && ((l != j) || (k != i)) && (diagonals || (l == j) || (k == i)))
					{
						if ((false == use_mask) || (mask_value == mask->GetPixel(k, l)))
						{
							int pix = pixeldata->GetPixel(k, l);
							if (pix != identifier)
							{
								if (pix > 0)
								{
									if (false == mode)
									{
										EdgePixels.insert(NewEdgePixels.begin(), NewEdgePixels.end());
										NewEdgePixels.clear();
										return (pix);
									}
								}
								else
								{
									surrounded = false;
									unsigned char gpix = gradientdata->GetPixel(k, l);
									if (gpix <= value)
									{
										ret = identifier;
										pixeldata->SetPixel(k, l, identifier);
										PointPair addpix;
										addpix.x = k;
										addpix.y = l;
										//										NewEdgePixels.insert(XY2Pos(addpix));
										NewEdgePixels.push_back(XY2Pos(addpix));
										size++;
										if (k < boundingbox.x0)
											boundingbox.x0 = k;
										if (k > boundingbox.x1)
											boundingbox.x1 = k;
										if (l < boundingbox.y0)
											boundingbox.y0 = l;
										if (l > boundingbox.y1)
											boundingbox.y1 = l;
									}
								}
							}
						}
					}
				}
			}
			if (false == surrounded)
			{
				PointPair addpix;
				addpix.x = i;
				addpix.y = j;
				//				NewEdgePixels.insert(XY2Pos(addpix));
				NewEdgePixels.push_back(XY2Pos(addpix));
			}
		}
		EdgePixels.clear();
		EdgePixels.insert(NewEdgePixels.begin(), NewEdgePixels.end());
		NewEdgePixels.clear();
	}
	if (0 == ret)
	{
		level_complete = value;
	}
	return ret;
}

int SuperPixel::GetPrevSize()
{
	return prevsize;
}

int SuperPixel::GetSize()
{
	return size;
}

int SuperPixel::GetSetSize()
{
	int ret = 0;
	SuperPixel* current = GetHead();
	while (NULL != current)
	{
		++ret;
		current = current->GetNext();
	}
	return ret;
}

bool SuperPixel::SetPrevSize()
{
	prevsize = size;
	return true;
}

bool SuperPixel::SetWindow(RectQuad w)
{
	boundingbox = w;
	return true;
}

bool SuperPixel::SetWorkspace(WorkSpace* ws)
{
	workspace = ws;
	return true;
}

int SuperPixel::GetIdentifier()
{
	return identifier;
}

SuperPixel* SuperPixel::GetNext()
{
	return next;
}

SPixelData* SuperPixel::GetPixelData()
{
	return pixeldata;
}

SuperPixel* SuperPixel::GetPrev()
{
	return prev;
}

bool SuperPixel::SetNext(SuperPixel* n)
{
	next = n;
	return true;
}

bool SuperPixel::SetPrev(SuperPixel* p)
{
	prev = p;
	return true;
}

SuperPixel* SuperPixel::GetHead()
{
	SuperPixel* head;
	if (NULL == prev)
	{
		head = this;
	}
	else {
		head = prev;
		while (NULL != head->GetPrev())
		{
			head = head->GetPrev();
		}
	}
	return head;
}

SuperPixel* SuperPixel::GetTail()
{
	SuperPixel* tail;
	if (NULL == next)
	{
		tail = this;
	}
	else {
		tail = next;
		while (NULL != tail->GetNext())
		{
			tail = tail->GetNext();
		}
	}
	return tail;
}

SuperPixelType SuperPixel::GetType()
{
	return type;
}

SuperPixel* SuperPixel::GetByIdentifier(int key)
{
	if (workspace != NULL)
	{
		return workspace->GetByIdentifier(key);
	}
	SuperPixel* current;
	current = GetHead();
	while (NULL != current)
	{
		if (current->GetIdentifier() == key)
		{
			return current;
		}
		current = current->GetNext();
	}
	return NULL;
}

bool SuperPixel::Absorb(SuperPixel* intersecting_pixel, bool remove, bool manage_neighbors)
{
	// "remove" indicates whether the intersecting_pixel should be deleted after being absorbed.
	//double channel_squared[3];
	Color other_color = intersecting_pixel->GetAveColor();
	int other_size = intersecting_pixel->GetSize();
	for (int chan = 0; chan < 3; ++chan)
	{
		//channel_squared[chan] = size * (AveColor.channel[chan] * AveColor.channel[chan]);

		//channel_squared[chan] += other_size * (other_color.channel[chan] * other_color.channel[chan]);

		AveColor.channel[chan] = sqrt((size * (AveColor.channel[chan] * AveColor.channel[chan]) + other_size * (other_color.channel[chan] * other_color.channel[chan])) / (size + other_size));
	}
	int int_id = intersecting_pixel->GetIdentifier();
	size += intersecting_pixel->GetSize();
	prevsize += intersecting_pixel->GetPrevSize();
	RectQuad wind = intersecting_pixel->GetWindow();
	if (wind.x0 < boundingbox.x0)
	{
		boundingbox.x0 = wind.x0;
	}
	if (wind.x1 > boundingbox.x1)
	{
		boundingbox.x1 = wind.x1;
	}
	if (wind.y0 < boundingbox.y0)
	{
		boundingbox.y0 = wind.y0;
	}
	if (wind.y1 > boundingbox.y1)
	{
		boundingbox.y1 = wind.y1;
	}
	for (int j = boundingbox.y0; j <= boundingbox.y1; j++)
	{
		for (int i = boundingbox.x0; i <= boundingbox.x1; i++)
		{
			if (int_id == pixeldata->GetPixel(i, j))
			{
				pixeldata->SetPixel(i, j, identifier);
			}
		}
	}
	level_complete = 0;
	// Update EdgePixels.
	FindEdgePixels();
	// EdgePixels.insert(intersecting_pixel->GetEdgePixels()->begin(), intersecting_pixel->GetEdgePixels()->end());
	if (manage_neighbors)
	{
		Neighbors.insert(intersecting_pixel->GetNeighbors()->begin(), intersecting_pixel->GetNeighbors()->end());
		for (std::set<int>::iterator iter = Neighbors.begin(); iter != Neighbors.end();)
		{
			if ((*iter == identifier) || (*iter == intersecting_pixel->GetIdentifier()))
			{
				iter = Neighbors.erase(iter);
			}
			else
			{
				SuperPixel* temp = GetByIdentifier(*iter);
				if (NULL != temp)
				{
					temp->RemoveNeighbor(intersecting_pixel->GetIdentifier());
					temp->AddNeighbor(identifier);
				}
				++iter;
			}
		}
	}
	if (remove)
	{
		delete intersecting_pixel;
	}
	else {
		intersecting_pixel->GetEdgePixels()->clear();
		intersecting_pixel->EdgePixelsCurrent = false;
	}
	return true;
}

RectQuad SuperPixel::GetWindow()
{
	return boundingbox;
}

WorkSpace* SuperPixel::GetWorkspace()
{
	return workspace;
}

PointPair SuperPixel::GetSeed()
{
	return seed;
}

bool SuperPixel::FindEdgePixels()
{
	bool left_surrounded = false;
	int height = pixeldata->GetHeight();  // *** Should these both just be image_height and image_width?
	int width = pixeldata->GetWidth();
	EdgePixels.clear();
	for (int j = boundingbox.y0; j <= boundingbox.y1; j++)
	{
		for (int i = boundingbox.x0; i <= boundingbox.x1; i++)
		{
			if (identifier == pixeldata->GetPixel(i, j))
			{
				bool edge;
				if ((i == 0) || (i == (width - 1)) || (j == 0) || (j == (height - 1)))  // Account for pixels on the edge of the image.
				{
					edge = true;
				}
				else {
					edge = false;
					for (int l = (j - 1); l <= (j + 1); l++)
					{
						if ((l >= 0) && (l < image_height))
						{
							if (left_surrounded)
							{
								int k = i + 1; // Can skip i-1 and i, since we know the pixel to the left (i-1) is surrounded.
								if (k < image_width)
								{
									if (false == edge)
									{
										if (identifier != pixeldata->GetPixel(k, l))
										{
											edge = true;
										}
									}
								}

							}
							else {
								for (int k = (i - 1); k <= (i + 1); k++)// Do full examination, because the pixel to the left is not surrounded.
								{
									if ((k >= 0) && (k < image_width))
									{
										if ((false == edge) && ((l != j) || (k != i)))
										{
											if (identifier != pixeldata->GetPixel(k, l))
											{
												edge = true;
											}
										}
									}
								}
							}
						}
					}
				}
				if (edge)
				{
					PointPair temp;
					temp.x = i;
					temp.y = j;
					EdgePixels.insert(XY2Pos(temp));
					// Look for whether this is a convergence point of three or more identifiers.  Point is to the upper left of this pixel, so examine pixels to left, above, and above left.
					bool diff_found = false;
					bool meeting_point = false;
					int first_diff = 0;
					for (int k = (i - 1); k <= i; k++)
					{
						for (int l = (j - 1); (l <= j) && ((k != i) || (l != j)); l++)
						{
							int local_id;
							if ((k >= 0) && (l >= 0))
							{
								local_id = pixeldata->GetPixel(k, l);
							}
							else {
								local_id = 0;
							}
							if (diff_found)
							{
								if ((local_id != identifier) && (local_id != first_diff))
								{
									meeting_point = true;
								}
							}
							else {
								if (local_id != identifier)
								{
									first_diff = local_id;
									diff_found = true;
								}
							}
						}
					}
					if (meeting_point)
					{
						pixeldata->SetMeetingPoint(temp);
					}
				}
				left_surrounded = !edge;
			}
			else {
				left_surrounded = false;
			}
		}
		left_surrounded = false;
	}
	EdgePixelsCurrent = true;
	return true;
}

bool SuperPixel::FindPaths(bool use_meeting_points, bool polygon, bool fine)
{
	Path* current = NULL;
	if (NULL != path_list_head)
	{
		current = path_list_head->GetHead();
	}
	while (NULL != current)
	{
		Path* next_path = current->GetNext();
		delete current;
		current = next_path;
	}
	path_list_head = NULL;
	path_list_tail = NULL;

	std::set<int>::iterator it;
	if (false == EdgePixelsCurrent)
	{
		FindEdgePixels();
	}
	if (EdgePixels.size() == 0)
	{
		return false;
	}

	/*
	it = EdgePixels.begin();
	PointPair start = Pos2XY(*it);
	Path* p = new Path(this, start, use_meeting_points, polygon, fine);
	if (NULL != p)  // *** p may have a zero value for PointSet.  Can't assume that the first point will return anything.  At some, point there may be none.  Fold this into remaining process below.
	{
		path_list_head = p;
		path_list_tail = p;
	}  // *** What if NULL == p?  Need to return false, right?
	std::set<int> p_set = p->GetPointSet();
	std::set<int> remaining;
	std::set_difference(EdgePixels.begin(), EdgePixels.end(), p_set.begin(), p_set.end(), std::inserter(remaining, remaining.begin()));

	*/

	Path* p = NULL;
	std::set<int> p_set;
	std::set<int> remaining;
	p_set.clear();
	remaining.clear();
	remaining.insert(EdgePixels.begin(), EdgePixels.end());

	int last_remaining_size = 0;
	while (remaining.size() != last_remaining_size)
	{
		last_remaining_size = remaining.size();
		for (it = remaining.begin(); it != remaining.end(); it++)
		{
			if ((NULL == path_list_head) || (NULL == path_list_head->GetByEdgePixel(Pos2XY(*it))))
			{
				p = new Path(this, Pos2XY(*it), use_meeting_points, polygon, fine);
				if (NULL == p)
				{
					throw std::runtime_error("No path found for edge pixel.\n");
					return false;
				}
				if (p->GetPointSet().size() > 3) 
				{
					if (NULL == path_list_head)
					{
						path_list_head = p;
						path_list_tail = p;
					}
					else {
						path_list_tail->InsertNext(p);
						path_list_tail = p;
					}
					p_set = p->GetPointSet();
					std::set<int> r_set(remaining);
					remaining.clear();
					std::set_difference(r_set.begin(), r_set.end(), p_set.begin(), p_set.end(), std::inserter(remaining, remaining.begin()));
					break;
				}
			}
		}
	}
	// Set internal paths to be sub-paths of outer paths.
	current = path_list_head;
	while (NULL != current)
	{
		Path* subsequent_path = current->GetNext();
		if (false == current->GetForward())
		{
			RectQuad sub_box = current->GetBox();
			p = path_list_head;
			while (NULL != p)
			{
				if (p->GetForward())
				{
					RectQuad forward_box = p->GetBox();
					if ((forward_box.x0 <= sub_box.x0) && (forward_box.x1 >= sub_box.x1) &&
						(forward_box.y0 <= sub_box.y0) && (forward_box.y1 >= sub_box.y1))
					{
						Path* c_prev = current->GetPrev();
						Path* c_next = subsequent_path;
						if (NULL != c_prev)
						{
							c_prev->SetNext(c_next);
						}
						if (NULL != c_next)
						{
							c_next->SetPrev(c_prev);
						}
						//current->SetNext(NULL);  // Don't need, moved to InsertSubpath.
						//current->SetPrev(NULL);
						p->InsertSubpath(current);
						path_list_head = p->GetHead();
						path_list_tail = p->GetTail();
						p = p->GetTail(); // Ensure that this won't loop again.
					}
				}
				p = p->GetNext();
			}
		}
		current = subsequent_path;
	}
	return true;
}

bool SuperPixel::GenerateContrastImage(SPixelData* pdata, int radius)
{
	if (radius < 0)
	{
		radius = 0;
	}
	ImageData* img = NULL;
	SuperPixel* current = NULL;
	std::set<Contrast_Histogram> histogram_index;
	std::set<int> present_identifiers;
	unsigned int* hist_entries;
	Contrast_Histogram current_hist;
	int* disc_chord = NULL;
	int wx, wy, img_x, img_y, disc_index, pix_value, local_identifier, location;
	double channel_squared[3] = { 0, 0, 0 };
	unsigned long long count = 0;
	std::set<Contrast_Histogram>::iterator it;
	Color local_color;
	RectQuad bbox;
	bbox.x0 = boundingbox.x0 - CONTRAST_BOX_MARGIN;
	bbox.y0 = boundingbox.y0 - CONTRAST_BOX_MARGIN;
	bbox.x1 = boundingbox.x1 + CONTRAST_BOX_MARGIN;
	bbox.y1 = boundingbox.y1 + CONTRAST_BOX_MARGIN;
	int width = 1 + (bbox.x1 - bbox.x0);
	int height = 1 + (bbox.y1 - bbox.y0);
	int img_width = pdata->GetWidth();
	int img_height = pdata->GetHeight();

	// Prepare the hist_entries and histogram_index, which hold the actual histogram and a set that translates a SuperPixel identifier into the location in the histogram.
	// Also prepare the present_identifiers, which is a set of those identifiers in the histogram with more than zero entries.
	hist_entries = (unsigned int*)malloc(sizeof(unsigned int) * Count());
	if (NULL == hist_entries)
	{
		throw(std::runtime_error("Unable to allocate memory for hist_entries for GenerateContrastImage.\n"));
	}

	histogram_index.clear();
	present_identifiers.clear();
	img = new ImageData(NULL, width, height, 3, true);
	current = GetHead();
	while (NULL != current)
	{
		current_hist.identifier = current->GetIdentifier();
		current_hist.location = count;
		current_hist.sp = current;
		histogram_index.insert(current_hist);
		hist_entries[count] = 0;
		count++;
		current = current->GetNext();
	}

	// Prepare disc_chord, contains the horizontal size of the disc at each vertical offset from the origin.
	disc_chord = (int*)malloc(sizeof(int) * (radius + 1));
	if (NULL == disc_chord)
	{
		throw (std::runtime_error("Unable to allocate memory for disc_chord for GenerateContrastImage.\n"));
	}
	for (int i = 0; i <= radius; i++)
	{
		disc_chord[i] = sqrt(radius * radius - i * i);
	}

	// Calculate histogram for first pixel, at location bbox.x0, bbox.y0.
	for (wy = -radius; wy <= radius; ++wy)  // wx and wy are offsets within bbox.
	{
		disc_index = abs(wy);
		img_y = wy + bbox.y0;
		if ((img_y >= 0) && (img_y < img_height))
		{
			for (wx = -disc_chord[disc_index]; wx <= disc_chord[disc_index]; ++wx)
			{
				img_x = wx + bbox.x0;
				if ((img_x >= 0) && (img_x < img_width))
				{
					pix_value = pdata->GetPixel(img_x, img_y);
					if (pix_value > 0)
					{
						current_hist.identifier = pix_value;
						it = histogram_index.find(current_hist);
						if (it != histogram_index.end())
						{
							current_hist = *it;
							hist_entries[current_hist.location]++;  // Increase histogram count for this SuperPixel.
							present_identifiers.insert(pix_value);  // Add this encountered SuperPixel to list of ones in the current histogram.
						}
					}
				}
			}
		}
	}
	// Calculate color for location bbox.x0, bbox.y0.
	//for (int chan = 0; chan < 3; ++chan)
	//{
	//	channel_squared[chan] = 0;
	//}
	//count = 0;
	//for (std::set<int>::iterator present_it = present_identifiers.begin(); present_it != present_identifiers.end(); ++present_it)  // Loop through the SuperPixels present here.
	//{
	//	local_identifier = *present_it;
	//	current_hist.identifier = local_identifier;
	//	it = histogram_index.find(current_hist);
	//	if (it != histogram_index.end())
	//	{
	//		current_hist = *it;
	//		int local_count = hist_entries[current_hist.location];
	//		local_color = current_hist.sp->GetAveColor();
	//		for (int chan = 0; chan < 3; ++chan)
	//		{
	//			channel_squared[chan] += local_count * (local_color.channel[chan] * local_color.channel[chan]);
	//		}
	//		count += local_count;
	//	}
	//}
	//if (count > 0)
	//{
	//	for (int chan = 0; chan < 3; ++chan)
	//	{
	//		local_color.channel[chan] = sqrt(channel_squared[chan] / count);
	//	}
	//	img->SetPixel(bbox.x0, bbox.y0, local_color);
	//}
	int direction = 1;  // 1 for moving to the right, -1 for moving to the left.
	int i = bbox.x0;
	int j = bbox.y0;
	while (j <= bbox.y1)
	{
		// Now do a horizontal strip.
		if (((direction > 0) && (i > bbox.x0)) || ((direction < 0) && (i < bbox.x1))) // If i == bbox.x0 or bbox.x1, then the histogram may already have been filled out, no update needed.
		{
			for (wy = -radius; wy <= radius; wy++)
			{
				img_y = j + wy;
				if ((img_y >= 0) && (img_y < img_height))
				{
					disc_index = abs(wy);
					img_x = (i + direction * disc_chord[disc_index]);  // Add to the histogram.  Newly in the disc.
					if ((img_x >= 0) && (img_x < img_width))
					{
						pix_value = pdata->GetPixel(img_x, img_y);
						if (pix_value > 0)
						{
							current_hist.identifier = pix_value;
							it = histogram_index.find(current_hist);
							if (it != histogram_index.end())
							{
								current_hist = *it;
								hist_entries[current_hist.location]++;  // Increase histogram count for this SuperPixel.
								present_identifiers.insert(pix_value);  // Add this encountered SuperPixel to list of ones in the current histogram.
							}
						}
					}
					img_x = (i - direction * disc_chord[disc_index] - direction);  // Subtract from the histogram.  Just to the outside of the disc.
					if ((img_x >= 0) && (img_x < img_width))
					{
						pix_value = pdata->GetPixel(img_x, img_y);
						if (pix_value > 0)
						{
							current_hist.identifier = pix_value;
							it = histogram_index.find(current_hist);
							if (it != histogram_index.end())
							{
								current_hist = *it;
								hist_entries[current_hist.location]--;  // Decrease histogram count for this SuperPixel.
								if (0 == hist_entries[current_hist.location])
								{
									present_identifiers.erase(pix_value);  // Remove this encountered SuperPixel from the list of ones in the current histogram.
								}
							}
						}
					}
				}
			}
		}
		// Calculate color for location i, j.
		for (int chan = 0; chan < 3; ++chan)
		{
			channel_squared[chan] = 0;
		}
		count = 0;
		for (std::set<int>::iterator present_it = present_identifiers.begin(); present_it != present_identifiers.end(); ++present_it)  // Loop through the SuperPixels present here.
		{
			local_identifier = *present_it;
			current_hist.identifier = local_identifier;
			it = histogram_index.find(current_hist);
			if (it != histogram_index.end())
			{
				current_hist = *it;
				int local_count = hist_entries[current_hist.location];
				local_color = current_hist.sp->GetAveColor();
				for (int chan = 0; chan < 3; ++chan)
				{
					channel_squared[chan] += local_count * (local_color.channel[chan] * local_color.channel[chan]);
				}
				count += local_count;
			}
		}
		if (count > 0)
		{
			for (int chan = 0; chan < 3; ++chan)
			{
				local_color.channel[chan] = sqrt(channel_squared[chan] / count);
			}
			local_color = CalculateContrastColor(local_color);
			img->SetPixel(i - bbox.x0, j - bbox.y0, local_color);
		}
		// Move horizontally
		i += direction;
		if (((i < bbox.x0) && (direction < 0)) || ((i > bbox.x1) && (direction > 0))) // Need to drop down and reverse directions.
		{
			j++;
			i = i - direction;
			direction = -direction;
			if (j <= bbox.y1)
			{
				for (wx = -radius; wx <= radius; ++wx)
				{
					img_x = wx + i;
					if ((img_x >= 0) && (img_x < img_width))
					{
						disc_index = abs(wx);
						img_y = (j + disc_chord[disc_index]);  // Add to the histogram.  Bottommost in the disc.
						if ((img_y >= 0) && (img_y < img_height))
						{
							pix_value = pdata->GetPixel(img_x, img_y);
							if (pix_value > 0)
							{
								current_hist.identifier = pix_value;
								it = histogram_index.find(current_hist);
								if (it != histogram_index.end())
								{
									current_hist = *it;
									hist_entries[current_hist.location]++;  // Increase histogram count for this SuperPixel.
									present_identifiers.insert(pix_value);  // Add this encountered SuperPixel to list of ones in the current histogram.
								}
							}
						}
						img_y = (j - disc_chord[disc_index] - 1);  // Subtract from the histogram.  Just above the topmost in the disc.
						if ((img_y >= 0) && (img_y < img_height))
						{
							pix_value = pdata->GetPixel(img_x, img_y);
							if (pix_value > 0)
							{
								current_hist.identifier = pix_value;
								it = histogram_index.find(current_hist);
								if (it != histogram_index.end())
								{
									current_hist = *it;
									hist_entries[current_hist.location]--;  // Decrease histogram count for this SuperPixel.
									if (0 == hist_entries[current_hist.location])
									{
										present_identifiers.erase(pix_value);  // Remove this encountered SuperPixel from the list of ones in the current histogram.
									}
								}
							}
						}
					}
				}
			}
		}
	}

	fill_image.clear();
	if (false == img->CollapseWideData(true))
	{
		throw std::runtime_error("Failed to collapse data_wide in GenerateContrastImage.\n");
	}
	if (0 == stbi_write_png_to_func(write_png_to_mem, this, width, height, 3, img->GetData(), 3 * width))
	{
		throw std::runtime_error("Unable to calculate fill_image.\n");
	}

	free(disc_chord);
	delete(img);
	free(hist_entries);
	return true;
}

Color SuperPixel::CalculateContrastColor(Color opposing)
{
	ImageData* img = workspace->GetImage();
	Color ret = img->CIELABconvert(AveColor);
	Color op = img->CIELABconvert(opposing);
	float post = ret.channel[0] - (op.channel[0] - ret.channel[0]) / 2;;
	if (post < 0)
	{
		post = 0;
	}
	else if (post > 100)
	{
		post = 100;
	}
	ret.channel[0] = post;
	ret = img->RGBconvert(ret);

	return ret;
}

SuperPixel* SuperPixel::DuplicateSuperPixelSet(SPixelData* sd) // Creates a duplicate set of SuperPixels, using the given SPixelData if it is provided.
{
	SuperPixel* head = NULL;
	SuperPixel* current = GetHead();
	SPixelData* data = sd;
	if (NULL == sd)
	{
		data = current->GetPixelData();
	}
	else {
		data = sd;
	}
	WorkSpace* wspace = current->GetWorkspace();
	SuperPixelType t = current->GetType();
	head = new SuperPixel(*current, current->GetGradient(), data, NULL, NULL, wspace, t);
	current = current->GetNext();
	SuperPixel* LocalPrev = head;
	while (NULL != current)
	{
		SuperPixel* LocalSP = new SuperPixel(*current, current->GetGradient(), data, NULL, LocalPrev, wspace, t);
		LocalPrev->SetNext(LocalSP);
		current = current->GetNext();
		LocalPrev = LocalSP;
	}
	return head;
}

bool SuperPixel::UpdateWorkspaceSuperPixelSet(WorkSpace* ws)
{
	SuperPixel* current = GetHead();
	while (NULL != current)
	{
		current->SetWorkspace(ws);
		current = current->GetNext();
	}
	return true;
}

bool SuperPixel::DeletePathList()
{
	Path* current_path = NULL;
	Path* next_path = path_list_head;
	while (NULL != next_path)
	{
		current_path = next_path;
		next_path = current_path->GetNext();
		delete current_path;
	}
	path_list_head = NULL;
	path_list_tail = NULL;
	EdgePixelsCurrent = false;
	EdgePixels.clear();
	return true;
}

bool SuperPixel::CalculateSize()
{
	size = 0;
	for (int j = boundingbox.y0; j <= boundingbox.y1; j++)
	{
		for (int i = boundingbox.x0; i <= boundingbox.x1; i++)
		{
			if (identifier == pixeldata->GetPixel(i, j))
			{
				size++;
			}
		}
	}
	return true;
}

bool SuperPixel::Reset()
{
	pixeldata->SetPixel(seed.x, seed.y, identifier);
	EdgePixels.clear();
	EdgePixels.insert(XY2Pos(seed));
	EdgePixelsCurrent = true;
	boundingbox.x0 = seed.x;
	boundingbox.y0 = seed.y;
	boundingbox.x1 = seed.x;
	boundingbox.y1 = seed.y;
	size = 1;
	prevsize = 1;
	level_complete = 0;
	AveError = 0;
	return true;
}

Color SuperPixel::GetAveColor()
{
	return AveColor;
}

float SuperPixel::GetAveError()
{
	return AveError;
}

std::set<int>* SuperPixel::GetEdgePixels()
{
	return &EdgePixels;
}

std::string SuperPixel::GetFillImage()
{
	return fill_image;
}

GradData* SuperPixel::GetGradient()
{
	return gradientdata;
}

PointPair SuperPixel::Pos2XY(int pos)
{
	PointPair ret;
	ret.y = (int)(pos / image_width);
	ret.x = pos - (ret.y * image_width);
	return ret;
}

int SuperPixel::XY2Pos(PointPair xy)
{
	int ret;
	ret = (xy.y * image_width + xy.x);
	return ret;
}

FloatPointPair SuperPixel::GetCentroid()
{
	return centroid;
}

int SuperPixel::GetColorBucket()
{
	return color_bucket;
}

float SuperPixel::ColorDifference(Color c1, Color c2)
{
	float distance;
	float diff_1, diff_2, diff_3;

	diff_1 = c1.channel[0] - c2.channel[0];
	diff_2 = c1.channel[1] - c2.channel[1];
	diff_3 = c1.channel[2] - c2.channel[2];
	distance = sqrt(diff_1 * diff_1 + diff_2 * diff_2 + diff_3 * diff_3);
	return distance;
}

PointPair SuperPixel::Split(ImageData* image, bool diagonals)
{
	PointPair ret;
	ret.x = 0;
	ret.y = 0;
	float diff = 0;
	unsigned char min = 255;
	int count = 1;
	SPixelData* seedpixeldata;
	SuperPixel* head = NULL;
	SuperPixel* tail = NULL;
	int dx, dy, nx, ny;
	Cell* cell = NULL;

	if ((boundingbox.x1 - boundingbox.x0) >= 70)
	{
		nx = 7;
		dx = (boundingbox.x1 - boundingbox.x0) / 7;
	}
	else {
		nx = (boundingbox.x1 - boundingbox.x0) / 10;
		if (nx < 1)
		{
			nx = 1;
		}
		dx = (boundingbox.x1 - boundingbox.x0) / nx;
	}
	if ((boundingbox.y1 - boundingbox.y0) >= 70)
	{
		ny = 7;
		dy = (boundingbox.y1 - boundingbox.y0) / 7;
	}
	else {
		ny = (boundingbox.y1 - boundingbox.y0) / 10;
		if (ny < 1)
		{
			ny = 1;
		}
		dy = (boundingbox.y1 - boundingbox.y0) / ny;
	}

	if ((nx < 2) && (ny < 2))
	{
		return ret;
	}

	seedpixeldata = new SPixelData(image_width, image_height);

	//First, set the SuperPixel with identifier '1' (through count variable) as the existing one.
	seedpixeldata->Reset();
	head = new SuperPixel(count, gradientdata, seedpixeldata, seed, NULL, NULL, workspace, SPType_Plain);
	tail = head;

	// Next, find potential sub-SuperPixels.
	for (int j = 0; j < ny; j++)
	{
		for (int i = 0; i < nx; i++)
		{
			int cx = boundingbox.x0 + i * dx + (dx / 2);
			int cy = boundingbox.y0 + j * dy + (dy / 2);
			cell = new Cell(gradientdata, seedpixeldata, cx, cy, (int)dx, (int)dy, 2);
			if (cell->FindSeed(pixeldata, identifier, diagonals)) // Use the existing identifier values in pixeldata as the active mask.
			{
				count++;
				PointPair prop_seed = cell->GetSeed();
				tail = new SuperPixel(count, gradientdata, seedpixeldata, prop_seed, NULL, tail, workspace, SPType_Plain);
			}
			delete cell;
		}
	}

	// Now, try out each new potential minimum.
	SuperPixel* current;
	current = head->GetNext();
	while (current != NULL)
	{
		seedpixeldata->Reset();
		head->Reset();
		current->Reset();
		for (int level = 0; level < 256; level++)
		{
			bool done = false;
			while (false == done)
			{
				done = true;
				int grow_ret = head->Grow(level, true, true, boundingbox, pixeldata, identifier, diagonals);
				if (grow_ret > 0)
				{
					done = false;
				}
				grow_ret = current->Grow(level, true, true, boundingbox, pixeldata, identifier, diagonals);
				if (grow_ret > 0)
				{
					done = false;
				}
			}
		}
		if ((head->GetSize() > (size / 10)) && (current->GetSize() > (size / 10)))
		{
			head->SetAveColor(image);
			current->SetAveColor(image);
			float inst_diff = ColorDifference(head->GetAveColor(), current->GetAveColor());
			if (inst_diff > diff)
			{
				diff = inst_diff;
				ret = current->GetSeed();
			}
		}
		current = current->GetNext();
	}
	current = head->GetNext();
	while (current != NULL)
	{
		delete current;
		current = head->GetNext();
	}
	delete head;
	delete seedpixeldata;
	return ret;
}


bool SuperPixel::SeparateDiscontinuousSuperPixel()
{
	// Step 1: See whether there are disconnected edge paths.
	if (NULL == path_list_head)
	{
		// Need to find paths.
		FindPaths(false, false, false);
	}
	if (path_list_head == path_list_tail)
	{
		// There is only one contiguous area.
		return true;
	}
	// Step 2:  Loop through any other paths that are forward paths (e.g. not interior paths).
	Path* current_path = path_list_head->GetNext();
	int point = 0;
	while (NULL != current_path)
	{
		Path* subsequent_path = current_path->GetNext();
		point = *(current_path->GetPointSet().begin());

		// Step 3: Generate new identifier.
		int new_id = 0;
		SuperPixel* current = GetHead();
		while (NULL != current)
		{
			if (current->GetIdentifier() > new_id)
			{
				new_id = current->GetIdentifier();
			}
			current = current->GetNext();
		}
		new_id++; // Next available number is the new identifier.

		// Step 4: Starting at point, flood-fill the existing pixeldata at point to be new identifier.
		if (false == pixeldata->FloodReplace(point, identifier, new_id))  //*** Here.  Flood fill goes beyond the bounding box.  Why? ***
		{
			throw std::runtime_error("Path point not in original set.\n");
			return false;
		}

		// Step 5: Calculate new bounding boxes.
		RectQuad orig_box = { 0,0,0,0 };
		RectQuad new_box = { 0,0,0,0 };
		bool orig_seen = false;
		bool new_seen = false;
		for (int j = boundingbox.y0; j <= boundingbox.y1; j++)
		{
			for (int i = boundingbox.x0; i <= boundingbox.x1; i++)
			{
				int value = pixeldata->GetPixel(i, j);
				if (identifier == value)
				{
					if (orig_seen)
					{
						if (i < orig_box.x0)
						{
							orig_box.x0 = i;
						}
						if (i > orig_box.x1)
						{
							orig_box.x1 = i;
						}
						if (j < orig_box.y0)
						{
							orig_box.y0 = j;
						}
						if (j > orig_box.y1)
						{
							orig_box.y1 = j;
						}
					}
					else {
						orig_box.x0 = i;
						orig_box.y0 = j;
						orig_box.x1 = i;
						orig_box.y1 = j;
						orig_seen = true;
					}
				}
				else if (new_id == value)
				{
					if (new_seen)
					{
						if (i < new_box.x0)
						{
							new_box.x0 = i;
						}
						if (i > new_box.x1)
						{
							new_box.x1 = i;
						}
						if (j < new_box.y0)
						{
							new_box.y0 = j;
						}
						if (j > new_box.y1)
						{
							new_box.y1 = j;
						}
					}
					else {
						new_box.x0 = i;
						new_box.y0 = j;
						new_box.x1 = i;
						new_box.y1 = j;
						new_seen = true;
					}
				}
			}
		}
		if (false == new_seen)
		{
			throw std::runtime_error("Failed to find new pixels in SeparateDiscontinuousSuperPixel.\n");
			return false;
		}

		// Step 6: Confirm that there are original identifier pixels left.
		if (false == orig_seen)
		{
			// Replace original entirely with the new one.
			SetIdentifierandBox(new_id, new_box);
		}
		else {
			// Step 7: Create new SuperPixel.
			PointPair seed;
			seed.y = point / pixeldata->GetWidth();
			seed.x = point - (seed.y * pixeldata->GetWidth());
			SuperPixel* newSP = NULL;
			if (NULL == gradientdata)
			{
				newSP = new SuperPixel(new_id, pixeldata, seed, NULL, GetTail(), workspace, type);
			}
			else {
				newSP = new SuperPixel(new_id, gradientdata, pixeldata, seed, NULL, GetTail(), workspace, type);
			}
			if (NULL == newSP)
			{
				throw std::runtime_error("Failed to create new SuperPixel in SeparateDiscontinuousSuperPixel.\n");
				return false;
			}
			newSP->SetWindow(new_box);
			boundingbox = orig_box;

			// Step 8: Determine which paths (with embedded sub-paths) go with the new SuperPixel.
			Path* local_path = GetPathHead();
			while (NULL != local_path)
			{
				Path* local_next_path = local_path->GetNext();
				std::set<int> points = local_path->GetPointSet();
				int pos = *points.begin();
				int y = pos / pixeldata->GetWidth();
				int x = pos - (y * pixeldata->GetWidth());
				int value = pixeldata->GetPixel(x, y);
				if (new_id == value) // This path needs to move over to the new SuperPixel.
				{
					if (local_path == subsequent_path)  // Not sure this can happen.
					{
						subsequent_path = subsequent_path->GetNext();
					}
					local_path->MoveSuperPixel(newSP);
				}
				local_path = local_next_path;
			}

			// Step 9: Set the color of the new SuperPixel to the same as the original.
			newSP->SetAveColor(AveColor);
		}
		current_path = subsequent_path;
	}
	return true;
}

SuperPixel* SuperPixel::Thin(GradData* graddat, SPixelData* pixdat, SuperPixel* prev, bool only_one, bool glitch3)
{
	bool changed = true;
	SuperPixel* ret = NULL;

	ret = new SuperPixel(*this, graddat, pixdat, NULL, prev, workspace, SPType_Skeleton);
	ret->FindEdgePixels();

	if (ret->EdgePixels.size() > 0)
	{
		while (changed)
		{
			changed = false;
			changed = ret->Thin_Subiteration(0, changed);
			changed = ret->Thin_Subiteration(1, changed);
		}
	}
	ret->FindSkeletonPoints();
	ret->skeleton_collection->FindSkeletonLinks();
	ret->skeleton_collection->RecalcDistances();
	ret->skeleton_collection->FindLongestPaths(only_one);
	size = EdgePixels.size();
	if (NULL != ret->path_list_head)
	{
		while (NULL != ret->path_list_head->GetTail())
		{
			delete(ret->path_list_head->GetTail());
		}
		delete(ret->path_list_head);
	}
	ret->path_list_head = NULL;
	ret->path_list_tail = NULL;
	bool begin_endpoint = false;
	bool end_endpoint = false;
	PointPair temp = { 0, 0 };
	std::set<std::vector<int>>::iterator path_iterator;
	std::set<std::vector<int>> paths = ret->skeleton_collection->GetLongestPaths();
	std::set<int> endpoints = ret->skeleton_collection->GetEndpoints();
	for (path_iterator = paths.begin(); path_iterator != paths.end(); ++path_iterator)
	{
		std::vector<int> long_path = *path_iterator;
		if (long_path.size() > 0)
		{
			if (endpoints.end() != endpoints.find(*long_path.begin()))
			{
				begin_endpoint = true;
			}
			else {
				begin_endpoint = false;
			}
			if (endpoints.end() != endpoints.find(long_path.back()))
			{
				end_endpoint = true;
			}
			else {
				end_endpoint = false;
			}
			if (NULL == ret->path_list_head)
			{
				ret->path_list_head = new Path(ret, temp, false, false, false, &long_path, begin_endpoint, end_endpoint, glitch3);
				ret->path_list_head->SetPrev(NULL);
				ret->path_list_head->SetNext(NULL);
				ret->path_list_tail = ret->path_list_head;
			}
			else {
				Path* prev_tail = ret->path_list_tail;
				ret->path_list_tail = new Path(ret, temp, false, false, false, &long_path, begin_endpoint, end_endpoint, glitch3);
				prev_tail->SetNext(ret->path_list_tail);
				ret->path_list_tail->SetPrev(prev_tail);
			}
		}
	}

	return ret;
}

bool SuperPixel::Thin_Subiteration(int n, bool global_changed)
{
	bool changed = false;
	std::set<int> remove_list;
	remove_list.clear();
	std::set<int>::iterator it;
	for (it = EdgePixels.begin(); it != EdgePixels.end(); ++it)
	{
		int pos = *it;
		PointPair xy = Pos2XY(pos);
		int i = xy.x;
		int j = xy.y;
		bool p2, p3, p4, p5, p6, p7, p8, p9;
		int B = 0;
		if (j > 0)  // p9, p2, p3
		{
			if (i > 0)
			{
				p9 = (identifier == pixeldata->GetPixel(i - 1, j - 1));
				if (p9)
				{
					B++;
				}
			}
			else {
				p9 = false;
			}
			if (i < (image_width - 1))
			{
				p3 = (identifier == pixeldata->GetPixel(i + 1, j - 1));
				if (p3)
				{
					B++;
				}
			}
			else {
				p3 = false;
			}
			p2 = (identifier == pixeldata->GetPixel(i, j - 1));
			if (p2)
			{
				B++;
			}
		}
		else {
			p9 = false;
			p2 = false;
			p3 = false;
		}

		if (j < (image_height - 1))  // p7, p6, p5
		{
			if (i > 0)
			{
				p7 = (identifier == pixeldata->GetPixel(i - 1, j + 1));
				if (p7)
				{
					B++;
				}
			}
			else {
				p7 = false;
			}
			if (i < (image_width - 1))
			{
				p5 = (identifier == pixeldata->GetPixel(i + 1, j + 1));
				if (p5)
				{
					B++;
				}
			}
			else {
				p5 = false;
			}
			p6 = (identifier == pixeldata->GetPixel(i, j + 1));
			if (p6)
			{
				B++;
			}
		}
		else {
			p5 = false;
			p6 = false;
			p7 = false;
		}

		if (i > 0)  // p8
		{
			p8 = (identifier == pixeldata->GetPixel(i - 1, j));
			if (p8)
			{
				B++;
			}
		}
		else {
			p8 = false;
		}

		if (i < (image_width - 1))  // p4
		{
			p4 = (identifier == pixeldata->GetPixel(i + 1, j));
			if (p4)
			{
				B++;
			}
		}
		else {
			p4 = false;
		}
		if ((B >= 2) && (B <= 6))
		{
			int A = 0;
			if (!p2 && p3)
			{
				A++;
			}
			if (!p3 && p4)
			{
				A++;
			}
			if (!p4 && p5)
			{
				A++;
			}
			if (!p5 && p6)
			{
				A++;
			}
			if (!p6 && p7)
			{
				A++;
			}
			if (!p7 && p8)
			{
				A++;
			}
			if (!p8 && p9)
			{
				A++;
			}
			if (!p9 && p2)
			{
				A++;
			}
			if (1 == A)
			{
				if (0 == n)
				{
					if (!p4 || !p6 || (!p2 && !p8))
					{
						remove_list.insert(pos);
						changed = true;
					}
				}
				else {
					if (!p2 || !p8 || (!p4 && !p6))
					{
						remove_list.insert(pos);
						changed = true;
					}
				}

			}
		}
	}
	if (changed)
	{
		for (it = remove_list.begin(); it != remove_list.end(); ++it)
		{
			int pos = *it;
			PointPair xy = Pos2XY(pos);
			int i = xy.x;
			int j = xy.y;
			EdgePixels.erase(pos);
			pixeldata->SetPixel(i, j, 0);
		}
		for (it = remove_list.begin(); it != remove_list.end(); ++it)
		{
			int pos = *it;
			PointPair xy = Pos2XY(pos);
			int i = xy.x;
			int j = xy.y;
			for (int wy = j - 1; wy <= (j + 1); wy++)
			{
				if ((wy >= 0) && (wy < image_height))
				{
					for (int wx = i - 1; wx <= (i + 1); wx++)
					{
						if ((wx >= 0) && (wx < image_width) && ((wx != i) || (wy != j))) // Exclude center point.
						{
							if (identifier == pixeldata->GetPixel(wx, wy))
							{
								PointPair new_xy;
								new_xy.x = wx;
								new_xy.y = wy;
								int new_pos = XY2Pos(new_xy);
								EdgePixels.insert(new_pos);
							}
						}
					}
				}
			}
		}
	}
	if ((0 == EdgePixels.size()) && (remove_list.size() > 0))  // In rare situations, all points are removed (when there is a 2x2 square).  In such case, set one point.
	{
		it = remove_list.begin();
		int pos = *it;
		PointPair xy = Pos2XY(pos);
		int i = xy.x;
		int j = xy.y;
		EdgePixels.insert(pos);
		pixeldata->SetPixel(i, j, identifier);
	}
	remove_list.clear();
	return (changed || global_changed);
}

bool SuperPixel::SetNeighbors()
{
	//  Future exploration: Maybe start with one border pixel at one edge of the bounding box, and walk the circuit.

	Neighbors.clear();
	FindEdgePixels();
	if (EdgePixels.size() > 0)
	{
		std::set<int>::iterator it;
		for (it = EdgePixels.begin(); it != EdgePixels.end(); ++it)
		{
			int pos = *it;
			PointPair xy = Pos2XY(pos);
			int i = xy.x;
			int j = xy.y;
			bool surrounded = true;
			for (int l = j - 1; l <= j + 1; l++)
			{
				for (int k = i - 1; k <= i + 1; k++)
				{
					if ((l >= 0) && (k >= 0) && (l < image_height) && (k < image_width) && ((l != j) || (k != i)))
					{
						int pix = pixeldata->GetPixel(k, l);
						if ((pix != identifier) && (pix > 0))
						{
							Neighbors.insert(pix);
						}
					}
				}
			}
		}
	}
	else {
		return false;
	}
	return true;
}

bool SuperPixel::SetPathHead(Path* path)
{
	path_list_head = path;
	return true;
}

bool SuperPixel::SetPathTail()
{
	path_list_tail = path_list_head->GetTail();
	return true;
}

std::set<int>* SuperPixel::GetNeighbors()
{
	return &Neighbors;
}


bool SuperPixel::CheckVertex(int v, bool check_neighbors)
{
	if (Vertices.count(v) > 0)
	{
		return true;
	}
	if (check_neighbors)
	{
		for (std::set<int>::iterator iter = Neighbors.begin(); iter != Neighbors.end(); ++iter)
		{
			SuperPixel* nsp = GetByIdentifier(*iter);
			if (NULL != nsp)
			{
				if (nsp->CheckVertex(v, false))
				{
					return true;
				}
			}
		}
	}
	return false;
}

std::vector<int> SuperPixel::write_data()
{
	int width = (boundingbox.x1 - boundingbox.x0 + 1);
	int height = (boundingbox.y1 - boundingbox.y0 + 1);
	int count;
	std::vector<int> ret;
	ret.clear();

	int i, j;
	for (j = 0; j < height; j++)
	{
		i = 0;
		count = 0;
		while ((i < width) && (identifier != pixeldata->GetPixel(i + boundingbox.x0, j + boundingbox.y0)))  // Always start with zero
		{
			count++;
			i++;
		}
		ret.push_back(count);
		while (i < width)
		{
			count = 0;
			while ((i < width) && (identifier == pixeldata->GetPixel(i + boundingbox.x0, j + boundingbox.y0)))  // Count ones
			{
				count++;
				i++;
			}
			if (count > 0)
			{
				ret.push_back(count);
			}
			else {
				break;
			}
			count = 0;
			while ((i < width) && (identifier != pixeldata->GetPixel(i + boundingbox.x0, j + boundingbox.y0)))  // Count zeros
			{
				count++;
				i++;
			}
			if (count > 0)
			{
				ret.push_back(count);
			}
			else {
				break;
			}
		}
	}
	return ret;
}

int SuperPixel::Count()
{
	int count = 0;
	SuperPixel* current;
	current = GetHead();
	while (NULL != current)
	{
		count++;
		current = current->GetNext();
	}
	return count;
}

bool SuperPixel::RemoveNeighbor(int ne)
{
	Neighbors.erase(ne);
	return true;
}

bool SuperPixel::FindSkeletonPoints()
{
	if (SPType_Skeleton != type)
	{
		throw std::runtime_error("Attempting to find SkeletonPoints in non-skeleton SuperPixel.\n");
		return false;
	}
	if (NULL != skeleton_collection)
	{
		delete skeleton_collection;
	}
	skeleton_collection = new SkeletonPointCollection(this);
	return skeleton_collection->FindSkeletonPoints();

	//int i, j, k, l, rot, transition_count;
	//PointPair p;
	//bool prev_state, local_bit;
	//unsigned char bit0, bit1, bit2;
	//ClearSkeletonPoints();
	//for (j = boundingbox.y0; j <= boundingbox.y1; j++)
	//{
	//	for (i = boundingbox.x0; i <= boundingbox.x1; i++)
	//	{
	//		if (identifier == pixeldata->GetPixel(i, j))
	//		{
	//			transition_count = 0;
	//			prev_state = true; // Start out with this as true, so first time checking rot=0 won't result in a false-to-true transition. 
	//			for (rot = 0; rot < 9; rot++)  // Rotate through top pixel clockwise back to the top pixel again.
	//			{
	//				bit0 = rot & 1;
	//				bit1 = (rot & 2) >> 1;
	//				bit2 = (rot & 4) >> 2;
	//				l = j + (bit2 * 2 - 1) * (-1) * ((bit1 & bit0) + bit1 - 1);  // sets l, relative to j, as -1 -1 0 1 1 1 0 -1.
	//				k = i + (bit2 * 2 - 1) * (-1) * (bit0 | bit1);  // sets k, relative to i, as 0 1 1 1 0 -1 -1 -1.
	//				local_bit = false;
	//				if ((k >= boundingbox.x0) && (k <= boundingbox.x1) && (l >= boundingbox.y0) && (l <= boundingbox.y1))
	//				{
	//					local_bit = (identifier == pixeldata->GetPixel(k, l));
	//				}
	//				if (local_bit && (false == prev_state)) // This is a false-to-true transition.
	//				{
	//					transition_count++;
	//				}
	//				prev_state = local_bit;
	//			}
	//			if (1 == transition_count) // Endpoint.
	//			{
	//				p.x = i;
	//				p.y = j;
	//				AddSkeletonEndpoint(p);
	//			}
	//			else if (transition_count > 2) // Intersection point.
	//			{
	//				p.x = i;
	//				p.y = j;
	//				AddSkeletonIntersection(p);
	//			}
	//		}
	//	}
	//}
	//return true;
}

bool SuperPixel::AddSkeletonEndpoint(PointPair p)
{
	if (SPType_Skeleton != type)
	{
		throw std::runtime_error("Attempting to add skeleton endpoint in non-skeleton SuperPixel.\n");
		return false;
	}
	skeleton_endpoint_set.insert(XY2Pos(p));
	return true;
}

bool SuperPixel::AddSkeletonIntersection(PointPair p)
{
	if (SPType_Skeleton != type)
	{
		throw std::runtime_error("Attempting to add skeleton intersection in non-skeleton SuperPixel.\n");
		return false;
	}
	skeleton_intersection_set.insert(XY2Pos(p));
	return true;
}

std::set<int> SuperPixel::GetSkeletonEndpoints()
{
	if (SPType_Skeleton != type)
	{
		throw std::runtime_error("Attempting to get skeleton endpoints in non-skeleton SuperPixel.\n");
	}
	return skeleton_collection->GetEndpoints();
}

std::set<int> SuperPixel::GetSkeletonIntersections()
{
	if (SPType_Skeleton != type)
	{
		throw std::runtime_error("Attempting to get skeleton intersections in non-skeleton SuperPixel.\n");
	}
	return skeleton_collection->GetIntersections();
}

bool SuperPixel::AddNeighbor(int ne)
{
	Neighbors.insert(ne);
	return true;
}

bool SuperPixel::AddVertex(int v)
{
	Vertices.insert(v);
	return true;
}

bool SuperPixel::ClearSkeletonPoints()
{
	skeleton_endpoint_set.clear();
	skeleton_intersection_set.clear();
	return true;
}

std::string SuperPixel::SaveState()
{
	std::string out;
	std::ostringstream s;
	s << identifier << ", " << size << ", " << seed.x << ", " << seed.y << ", " << centroid.x << ", " << centroid.y << ", ";
	s << boundingbox.x0 << ", " << boundingbox.y0 << ", " << boundingbox.x1 << ", " << boundingbox.y1 << ", ";
	s << (int)AveColor.channel[0] << ", " << (int)AveColor.channel[1] << ", " << (int)AveColor.channel[2] << ", " << AveError << "\n";
	out = s.str();
	return out;
}

Color SuperPixel::SetAveColor(ImageData* image)
{
	double channel_squared[3] = { 0, 0, 0 };
	unsigned long long count = 0;
	unsigned long long cent_x = 0;
	unsigned long long cent_y = 0;

	// Calculate AveColor
	for (int j = boundingbox.y0; j <= boundingbox.y1; j++)
	{
		for (int i = boundingbox.x0; i <= boundingbox.x1; i++)
		{
			if (identifier == pixeldata->GetPixel(i, j))
			{
				count++;
				Color c = image->GetPixel(i, j);
				for (int chan = 0; chan < 3; ++chan)
				{
					channel_squared[chan] += c.channel[chan] * c.channel[chan];
				}
				cent_x += i;
				cent_y += j;
			}
		}
	}
	for (int chan = 0; chan < 3; ++chan)
	{
		AveColor.channel[chan] = sqrt(channel_squared[chan] / count);
	}
	centroid.x = (double)(cent_x) / count;
	centroid.y = (double)(cent_y) / count;

	// Calculate AveError
	AveError = 0;
	for (int j = boundingbox.y0; j <= boundingbox.y1; j++)
	{
		for (int i = boundingbox.x0; i <= boundingbox.x1; i++)
		{
			if (identifier == pixeldata->GetPixel(i, j))
			{
				Color c2 = image->GetPixel(i, j);
				AveError += ColorDifference(AveColor, c2);
			}
		}
	}
	AveError = AveError / (float)(count);
	return AveColor;
}

Color SuperPixel::SetAveColor(Color c)
{
	AveColor = c;
	return AveColor;
}

bool SuperPixel::SetColorBucket(int b)
{
	color_bucket = b;
	return true;
}

bool SuperPixel::SetFillImage(std::string fill)
{
	fill_image = fill;
	return true;
}

bool SuperPixel::SetSize(int s)
{
	size = s;
	return true;
}

bool SuperPixel::SetIdentifierandBox(int id, RectQuad box)
{
	identifier = id;
	boundingbox = box;
	return true;
}

Path* SuperPixel::GetPathHead()
{
	return path_list_head;
}

Path* SuperPixel::GetPathTail()
{
	return path_list_tail;
}

