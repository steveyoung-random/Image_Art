// Copyright (c) 2023-2025 Steve Young
// Licensed under the MIT License

#include "Path.h"
#include "Workspace.h"

Path::Path(SuperPixel* sp, PointPair start, bool use_meeting_points, bool create_polygons, bool fine, std::vector<int>* line, bool beg_ep, bool end_ep, bool glitch3)
{
	// closed indicates whether this is a closed polygon or an unclosed line.

	next = NULL;
	prev = NULL;
	subpath = NULL;
	superpixel = sp;
	begin_endpoint = beg_ep;
	end_endpoint = end_ep;
	potential_segs.clear();
	point_set.clear();
	points.clear();
	fine_points.clear();
	optimal_segments.clear();
	optimal_curve.clear();
	if (NULL == line)
	{
		closed = true;
	}
	else {
		closed = false;
	}

	Color sp_color = sp->GetAveColor();
	std::stringstream sstream;
	if (sp_color.channel[0] < 16)
	{
		sstream << "0";
	}
	sstream << std::hex << (int)sp_color.channel[0];
	if (sp_color.channel[1] < 16)
	{
		sstream << "0";
	}
	sstream << std::hex << (int)sp_color.channel[1];
	if (sp_color.channel[2] < 16)
	{
		sstream << "0";
	}
	sstream << std::hex << (int)sp_color.channel[2];
	color_string = sstream.str();

	boundingbox.x0 = 0;
	boundingbox.y0 = 0;
	boundingbox.x1 = 0;
	boundingbox.y1 = 0;

	if (closed)
	{
		CalcPoints(start, use_meeting_points);
		if (create_polygons)
		{
			CalcForwardSegments(use_meeting_points, fine);
			if (CalcOptimalSegments())
			{
				Curve(false);
			}
		}
	}
	else {
		std::vector<int>::iterator it;
		if (line->size() > 0)
		{
			PointPair p = superpixel->Pos2XY(*line->begin());
			boundingbox.x0 = p.x;
			boundingbox.y0 = p.y;
			boundingbox.x1 = p.x;
			boundingbox.y1 = p.y;
		}
		if (line->size() == 1)
		{
			SPixelData* Orig_SPixDat = superpixel->GetWorkspace()->GetPixeldata();
			if (NULL != Orig_SPixDat)
			{
				PointPair p = superpixel->Pos2XY(*line->begin());
				line->clear();
				FloatPointPair p0;
				p0.x = p.x;
				p0.y = p.y;
				FloatPointPair p1;
				p1.x = p0.x + 1;
				p1.y = p0.y;
				FloatPointPair axial_edge;
				axial_edge = Orig_SPixDat->AxialExtent(p0, p1, superpixel->GetIdentifier());
				p.x = axial_edge.x;
				p.y = axial_edge.y;
				boundingbox.x0 = p.x;
				boundingbox.y0 = p.y;
				float dist = abs(p0.x - p.x);
				if (dist >= 25)
				{
					p.x = p0.x - dist + 10;
					line->push_back(superpixel->XY2Pos(p));
				}
				else if (dist >= 10)
				{
					p.x = p0.x - dist + 5;
					line->push_back(superpixel->XY2Pos(p));
				}
				else if (dist > 1)
				{
					p.x = p0.x - dist + 1;
					line->push_back(superpixel->XY2Pos(p));
				}
				p.x = p0.x;
				line->push_back(superpixel->XY2Pos(p));

				p1.x = p0.x - 1;
				axial_edge = Orig_SPixDat->AxialExtent(p0, p1, superpixel->GetIdentifier());
				p.x = axial_edge.x;
				p.y = axial_edge.y;
				dist = abs(p0.x - p.x);
				if (dist >= 25)
				{
					p.x = p0.x + dist - 10;
					line->push_back(superpixel->XY2Pos(p));
				}
				else if (dist >= 10)
				{
					p.x = p0.x + dist - 5;
					line->push_back(superpixel->XY2Pos(p));
				}
				else if (dist > 1)
				{
					p.x = p0.x + dist - 1;
					line->push_back(superpixel->XY2Pos(p));
				}
				p.x = p0.x + dist;
				boundingbox.x1 = p.x;
				boundingbox.y1 = p.y;
			}
		}
		for (it = line->begin(); it != line->end(); ++it)
		{
			point_set.insert(*it);
			PointPair current = sp->Pos2XY(*it);
			points.push_back(current);
			fine_points.push_back(current);
			if (current.x < boundingbox.x0)
			{
				boundingbox.x0 = current.x;
			}
			if (current.x > boundingbox.x1)
			{
				boundingbox.x1 = current.x;
			}
			if (current.y < boundingbox.y0)
			{
				boundingbox.y0 = current.y;
			}
			if (current.y > boundingbox.y1)
			{
				boundingbox.y1 = current.y;
			}
		}
		CalcForwardSegments(false, false);
		if (CalcOptimalSegments())
		{
			Curve(glitch3);
		}

	}
	return;
}

Path::~Path()
{
	if (NULL != prev)
	{
		prev->SetNext(next);
	}
	if (NULL != next)
	{
		next->SetPrev(prev);
	}
	Path* current_subpath;
	if (NULL != subpath)
	{
		current_subpath = subpath->GetHead();
		while (NULL != current_subpath)
		{
			Path* temp = current_subpath->GetNext();
			delete current_subpath;
			current_subpath = temp;
		}
	}
}

std::vector<PointPair> Path::GetPoints()
{
	return points;
}

std::set<int> Path::GetPointSet()
{
	return point_set;
}

std::vector<PointPair> Path::GetFinePoints()
{
	return fine_points;
}

std::vector<PolygonVertex> Path::GetOptSegments()
{
	return optimal_segments;
}

std::vector<Corner> Path::GetCurve()
{
	return optimal_curve;
}

FloatPointPair Path::AdjustVertex(int prev_point, int current_point, int next_point)
{
	FloatPointPair ret;
	FloatPointPair alphabeta;
	FloatPointPair vertex;

	if ((prev_point == current_point) || (current_point == next_point)) // Two points are at the same place.
	{
		ret.x = path_points[current_point].x;
		ret.y = path_points[current_point].y;
		return ret;
	}

	// Lines are y1 = a1 + b1*x1, and y2 = a1 + b2*x2
	float a1, a2, b1, b2;
	alphabeta = LinearRegression(prev_point, current_point);
	a1 = alphabeta.x;
	b1 = alphabeta.y;
	alphabeta = LinearRegression(current_point, next_point);
	a2 = alphabeta.x;
	b2 = alphabeta.y;

	if ((a1 == 0) && (b1 == 0)) // Line 1 is vertical, or horizontal at y=0
	{
		if (abs(path_points[current_point].y) < EFFECTIVE_ZERO)
		{
			ret.x = path_points[current_point].x;
			ret.y = 0;
		}
		else {
			ret.x = path_points[current_point].x;
			ret.y = a2 + b2 * ret.x;
		}
	}
	else if ((a2 == 0) && (b2 == 0)) // Line 2 is vertical, or horizontal at y=0
	{
		if (abs(path_points[current_point].y) < EFFECTIVE_ZERO)
		{
			ret.x = path_points[current_point].x;
			ret.y = 0;
		}
		else {
			ret.x = path_points[current_point].x;
			ret.y = a1 + b1 * ret.x;
		}
	}
	else if (abs(b2 - b1) < EFFECTIVE_ZERO) // Parallel lines, could happen so account for it.
	{
		ret.x = path_points[current_point].x;
		ret.y = a1 + b1 * ret.x;
	}
	else { // General case.
		ret.x = (a2 - a1) / (b1 - b2);
		ret.y = a1 + b1 * ret.x;
	}
	vertex.x = path_points[current_point].x;
	vertex.y = path_points[current_point].y;
	if ((abs(vertex.x - ret.x) > 0.5) || (abs(vertex.y - ret.y) > 0.5))
	{
		ret = SquareNeighborhood(vertex, ret);
	}
	return ret;
}

RectQuad Path::GetBox()
{
	return boundingbox;
}

bool Path::GetForward()
{
	return forward;
}

bool Path::GetClosed()
{
	return closed;
}

std::string Path::GetColor()
{
	return color_string;
}

Path* Path::GetByEdgePixel(PointPair point)
{
	Path* ret = NULL;
	if ((point.x >= (boundingbox.x0 - 1)) && (point.x <= (boundingbox.x1 + 1)) && (point.y >= (boundingbox.y0 - 1)) && (point.y <= (boundingbox.y1 + 1)))
	{
		std::vector<PointPair>::iterator p_it;
		p_it = points.begin();
		for (p_it = points.begin(); p_it != points.end(); ++p_it)
		{
			PointPair p = *p_it;
			if ((p.x == point.x) && (p.y == point.y))
			{
				ret = this;
				return ret;
			}
			else if ((p.y == point.y) && (1 == abs(p.x - point.x)))
			{
				ret = this;
				return ret;
			}
			else if ((p.x == point.x) && (1 == abs(p.y - point.y)))
			{
				ret = this;
				return ret;
			}
		}
	}
	if ((NULL == ret) && (NULL != next))
	{
		ret = next->GetByEdgePixel(point);
	}
	return ret;
}

Path* Path::GetHead()
{
	Path* head;
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

Path* Path::GetNext()
{
	return next;
}

Path* Path::GetPrev()
{
	return prev;
}

Path* Path::GetTail()
{
	Path* tail;
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

Path* Path::GetSubpathHead()
{
	if (NULL != subpath)
	{
		return subpath->GetHead();
	}
	return NULL;
}

Path* Path::GetSubpathTail()
{
	if (NULL != subpath)
	{
		return subpath->GetTail();
	}
	return nullptr;
}

FloatPointPair Path::LinearRegression(int point1, int point2)
{
	FloatPointPair ret;  // Not used as a cartesian coordinate point, but instead x = alpha and y = beta for equation of a line.
	float sum_xy = 0;
	float sum_x2 = 0;
	float ave_x = 0;
	float ave_y = 0;
	int n = path_points.size();
	int count = 0;

	int i = point1;
	bool cont = true;
	while (cont) // First loop to calculate averages.
	{
		ave_x += path_points[i].x;
		ave_y += path_points[i].y;
		count++;
		if (i == point2)
		{
			cont = false;
		}
		i++;
		if (i >= n)
		{
			i = 0;
		}
	}
	ave_x = ave_x / count;
	ave_y = ave_y / count;

	i = point1;
	cont = true;
	while (cont) // Second loop to calculate other values.
	{
		sum_xy += (path_points[i].x - ave_x) * (path_points[i].y - ave_y);
		sum_x2 += ((path_points[i].x - ave_x) * (path_points[i].x - ave_x));
		if (i == point2)
		{
			cont = false;
		}
		i++;
		if (i >= n)
		{
			i = 0;
		}
	}
	if (sum_x2 < EFFECTIVE_ZERO)
	{
		ret.x = 0;
		ret.y = 0;
		return ret;
	}
	ret.y = sum_xy / sum_x2; // Beta
	ret.x = ave_y - ret.y * ave_x;  // Alpha

	return ret;
}

FloatPointPair Path::SquareNeighborhood(FloatPointPair vertex, FloatPointPair point)
{
	FloatPointPair ret;
	float alpha, beta;

	if (abs(point.x - vertex.x) < EFFECTIVE_ZERO) // Vertical
	{
		ret.x = vertex.x;
		if (point.y < vertex.y)
		{
			ret.y = vertex.y - 0.5;
		}
		else {
			ret.y = vertex.y + 0.5;
		}
		return ret;
	}

	beta = (point.y - vertex.y) / (point.x - vertex.x);
	alpha = point.y - beta * point.x;

	if (point.x < (vertex.x - 0.5))
	{
		if (beta > 1.0)  // Will intersect at top.
		{
			ret.y = vertex.y - 0.5;
			ret.x = (ret.y - alpha) / beta;
		}
		else if (beta < -1.0)  // Will intersect at bottom.
		{
			ret.y = vertex.y + 0.5;
			ret.x = (ret.y - alpha) / beta;
		}
		else {  // Will intersect on left side.
			ret.x = vertex.x - 0.5;
			ret.y = alpha + beta * ret.x;
		}
		return ret;
	}
	else if (point.x > (vertex.x + 0.5))
	{
		if (beta > 1.0)  // Will intersect at bottom.
		{
			ret.y = vertex.y + 0.5;
			ret.x = (ret.y - alpha) / beta;
		}
		else if (beta < -1.0)  // Will intersect at top.
		{
			ret.y = vertex.y - 0.5;
			ret.x = (ret.y - alpha) / beta;
		}
		else {  // Will intersect on right side.
			ret.x = vertex.x + 0.5;
			ret.y = alpha + beta * ret.x;
		}
		return ret;
	}

	// X not beyond square, so Y must be.
	if (point.y < vertex.y)
	{
		ret.y = vertex.y - 0.5;
	}
	else {
		ret.y = vertex.y + 0.5;
	}
	ret.x = (vertex.y - alpha) / beta;
	return ret;
}

bool Path::CalcPoints(PointPair start, bool use_meeting_points)
{
	std::set<int>* EdgePixels = superpixel->GetEdgePixels();
	//	GradData* gradient = superpixel->GetGradient();
	SPixelData* pixeldata = superpixel->GetPixelData();

	int identifier = superpixel->GetIdentifier();
	int width = pixeldata->GetWidth();
	int height = pixeldata->GetHeight();
	int neighbors = 0;
	int direction = 0;
	int initial_direction = 0;
	bool first_flag = true;
	int cumulative_directions = 0;
	point_set.clear();
	PointPair initial_point;
	PointPair offset;
	// Directions:
	// 0 - undefined
	// 1 - right
	// 2 - down
	// 3 - left
	// 4 - up

	boundingbox.x0 = 0;
	boundingbox.y0 = 0;
	boundingbox.x1 = 0;
	boundingbox.y1 = 0;

	PointPair current = start;
	bool ready = false;
	while (false == ready)
	{
		neighbors = pixeldata->GetNeighbors(identifier, current);
		if ((-1 == neighbors) || (0 == (neighbors & 85)))
		{
			// Valid result if only one pixel in this object.
			points.clear();
			//fine_points.clear();
			points.push_back(current);
			boundingbox.x0 = current.x;
			boundingbox.y0 = current.y;
			boundingbox.x1 = current.x;
			boundingbox.y1 = current.y;
			point_set.insert(current.y * width + current.x);
			//fine_points.push_back(current);
			//current.x++;
			//fine_points.push_back(current);
			//current.y++;
			//fine_points.push_back(current);
			//current.x--;
			//fine_points.push_back(current);
			return true;
		}
		if ((neighbors & 255) == 255)
		{
			throw std::runtime_error("Start pixel in pathfinding is entirely enclosed.\n");
			return false;
		}
		if ((neighbors & 85) == 85)  // All adjacent pixels occupied.
		{
			// No adjacent open pixels, so look for open diagonals, and move to a pixel that isn't in a corner.
			if (0 == (neighbors & 2))
			{
				current.y -= 1;
			}
			else if (0 == (neighbors & 8))
			{
				current.x += 1;
			}
			else if (0 == (neighbors & 32))
			{
				current.y += 1;
			}
			else {
				current.x -= 1;
			}
		}
		else
		{
			ready = true;
		}
	}
	points.clear();
	fine_points.clear();

	if (0 == (neighbors & 1)) // Open above
	{
		direction = 1;
	}
	else if (0 == (neighbors & 4))
	{
		direction = 2;
	}
	else if (0 == (neighbors & 16))
	{
		direction = 3;
	}
	else if (0 == (neighbors & 64))
	{
		direction = 4;
	}
	else {
		std::runtime_error("Logic problem in Path function.\n");
		return false;
	}

	// Begin populating the points collection.
	points.push_back(current);
	boundingbox.x0 = current.x;
	boundingbox.y0 = current.y;
	boundingbox.x1 = current.x;
	boundingbox.y1 = current.y;
	point_set.insert(current.y * width + current.x);
	initial_direction = direction;
	initial_point = current;
	// The fine points are calculated by moving along the outer edges of the pixels, always keeping the filled pixels to the right and open to the left.
	if (direction > 2)     // d1 +.5, d2 +.5, d3 -.5, d4 -.5  (all values for offset are +.5 to those shown to the left, so only integers are used).
	{
		offset.x = current.x;  // Offset represents point to the upper left of offset pixel value.
	}
	else {
		offset.x = current.x + 1;
	}
	if (2 == (direction & 2))  // Directions 2 and 3 (down and to the left).
	{
		offset.y = current.y + 1; // d1 -.5, d2 +.5, d3 +.5, d4 -.5 (all values for offset are +.5 to those shown to the left, so only integers are used).
	}
	else {
		offset.y = current.y;
	}
	fine_points.push_back(offset);

	while (first_flag || (current.x != initial_point.x) || (current.y != initial_point.y) || (direction != initial_direction))  // Main loop.
	{
		first_flag = false;
		neighbors = pixeldata->GetNeighbors(identifier, current);

		// Do some preliminary calculations for the direction-based rules below.
		int ahead, ahead_right;
		ahead = 1 << (2 * direction - 1);
		ahead_right = ahead * 2;
		if (ahead_right > 128)
		{
			ahead_right = 1;
		}
		int x_left, x_straight;
		int y_left, y_straight;

		x_straight = (direction & 1);			// d1 +1, d2 0, d3 -1 d4 0
		if (2 == (direction & 2))
		{
			x_straight = -x_straight;
		}
		if (direction > 2) 		             	// d1 +1, d2 +1, d3 -1, d4 -1
		{
			x_left = -1;
		}
		else {
			x_left = 1;
		}
		y_straight = x_left - x_straight;      // d1 0, d2 +1, d3 0, d4 -1
		y_left = y_straight - x_straight;      // d1 -1, d2 +1, d3 +1, d4 -1

		if (0 == (neighbors & ahead_right)) // Right turn
		{
			direction += 1;
			if (direction > 4)
			{
				direction = 1;
			}
			cumulative_directions++;
		}
		else if (0 == (neighbors & ahead)) // Straight
		{
			current.x = current.x + x_straight;
			current.y = current.y + y_straight;
			points.push_back(current);
			point_set.insert(current.y * width + current.x);
		}
		else  // Left turn
		{
			direction -= 1;
			if (direction < 1)
			{
				direction = 4;
			}
			cumulative_directions--;
			current.x = current.x + x_left;
			current.y = current.y + y_left;
			points.push_back(current);
			point_set.insert(current.y * width + current.x);
		}
		if (direction > 2)     // d1 +.5, d2 +.5, d3 -.5, d4 -.5 (all values for offset are +.5 to those shown to the left, so only integers are used).
		{
			offset.x = current.x;
		}
		else {
			offset.x = current.x + 1;
		}
		if (2 == (direction & 2))
		{
			offset.y = current.y + 1; // d1 -.5, d2 +.5, d3 +.5, d4 -.5 (all values for offset are +.5 to those shown to the left, so only integers are used).
		}
		else {
			offset.y = current.y;
		}
		fine_points.push_back(offset);

		if (current.x < boundingbox.x0)
		{
			boundingbox.x0 = current.x;
		}
		else if (current.x > boundingbox.x1)
		{
			boundingbox.x1 = current.x;
		}
		if (current.y < boundingbox.y0)
		{
			boundingbox.y0 = current.y;
		}
		else if (current.y > boundingbox.y1)
		{
			boundingbox.y1 = current.y;
		}
	}
	forward = (cumulative_directions > 0);
	if (fine_points.size() > 1)
	{
		PointPair p0 = fine_points[0];
		PointPair plast = fine_points[fine_points.size() - 1];
		if ((p0.x == plast.x) && (p0.y == plast.y))
		{
			fine_points.pop_back();
		}
	}
	return true;
}

bool Path::CalcForwardSegments(bool use_meeting_points, bool fine)
{
	int start = 0;
	if (false == closed)
	{
		use_meeting_points = false;
		fine = false;
	}
	if (fine)
	{
		path_points = fine_points;
	}
	else {
		path_points = points;
	}
	int n = path_points.size();
	// If meeting_points are being used, we need to make sure we start at a meeting_point, if there is one.
	if (use_meeting_points && fine)
	{
		int i = 0;
		while ((false == superpixel->GetPixelData()->GetMeetingPoint(path_points[i])) && (false == superpixel->CheckVertex(superpixel->XY2Pos(path_points[i]), true)))
		{
			i++;
			if (i >= n)
			{
				i = 0;
				break;
			}
		}
		start = i;
	}

	int direction; // 1 = right, 2 = down, 3 = left, 4 = up
	bool dir_count[4];
	potential_segs.clear();

	for (int counter = 0; counter < n; counter++)
	{
		int i = (counter + start) % n;

		Segments segments;
		segments.vertex = i;
		if (use_meeting_points && superpixel->GetPixelData()->GetMeetingPoint(path_points[i]))
		{
			segments.meeting_point = true;
		}
		else {
			segments.meeting_point = false;
		}
		segments.s.clear();
		potential_segs.push_back(segments);

		PointPair p0 = path_points[i];
		int j;
		if (closed)
		{
			j = (i + 1) % n;
		}
		else {
			j = i + 1;
		}
		for (direction = 0; direction < 4; direction++)
		{
			dir_count[direction] = false;
		}
		bool straight = true;
		if (j >= n)
		{
			straight = false;
		}
		int NonStraightCounter = 0; // Even if a few j endpoints aren't straight, further ones may be.
		while (straight) // As long as the segments are still straight, keep extending out j.
		{
			if (i != j)
			{
				int j1 = j - 1;
				if (j1 < 0) // To account for wrap-around.
				{
					j1 = n - 1;
				}
				PointPair p1 = path_points[j1]; // Initially, find direction for last segment before j.
				PointPair p2 = path_points[j];
				if (p1.y == p2.y) // Horizontal
				{
					if (p2.x > p1.x)
					{
						direction = 1;
					}
					else
					{
						direction = 3;
					}
				}
				else { // Vertical
					if (p2.y > p1.y)
					{
						direction = 2;
					}
					else {
						direction = 4;
					}
				}
				dir_count[direction - 1] = true; // Update direction count
				if (dir_count[0] && dir_count[1] && dir_count[2] && dir_count[3]) // If all directions are seen, not a possible segment.
				{
					straight = false;
				}
				else {
					float dx = p2.x - p0.x;
					float dy = p2.y - p0.y;
					if ((dx != 0) || (dy != 0))
					{
						float inv_d02 = 1.0 / sqrt(dx * dx + dy * dy);
						float penalty = 0;
						bool LocalNonStraight = false;
						for (int k = (i + 1) % n; straight && (k != j); k = (k + 1) % n) // Cover all intermediary points between i and j.  Yes, I know that k doesn't come between i and j in the alphabet.
						{
							p1 = path_points[k];  // p1 is a point between p0 and p2, but need to confirm it is close enough to the line between p0 and p2.
							float dist = abs(dx * (p1.y - p0.y) - dy * (p1.x - p0.x)) * inv_d02;  // Note, no check is done for p1 not being beyond the ends of the line segment.  Could introduce an inaccuracy.
							if (dist > LINE_DISTANCE)
							{
								NonStraightCounter++;
								LocalNonStraight = true;
								if (NonStraightCounter > 20)
								{
									straight = false;
								}
							}
							if (straight && (LocalNonStraight == false))
							{
								penalty += dist;  // This, or squared?
							}
						}
						if (straight)
						{
							if (LocalNonStraight == false) // Only true if none of the intermediate points are too far away from the line.
							{
								Segment seg;
								seg.j = j;
								seg.penalty = penalty;
								potential_segs[counter].s.push_back(seg);
								NonStraightCounter = 0; // Reset counter if a straight line is encountered.
							}
							if (use_meeting_points && fine)  // If j is a meeting_point or neighboring vertex, then we cannot allow a segment to jump over it, so it has to be the last segment with this i.
							{
								if ((superpixel->GetPixelData()->GetMeetingPoint(path_points[j])) || (superpixel->CheckVertex(superpixel->XY2Pos(path_points[j]), true)))
								{
									straight = false;  // Not bailing because the segment beyond isn't straight, just using this as a flag to keep the segment from going past this point.
								}
							}
						}
					}
					else { // i and j point to the same location.  Rare, but possible (maybe).
						straight = false;
					}
				}
			}
			else {
				straight = false;
			}
			if (closed)
			{
				j = (j + 1) % n;
			}
			else {
				j = j + 1;
				if (j >= n)
				{
					straight = false;
				}
			}
		}
	}
	return true;
}

bool Path::CalcOptimalSegments()
{
	if (false == CalcSegmentPenalties())
	{
		return false;
	}

	int start = potential_segs[0].vertex;  // Offset for first vertex.
	int n = potential_segs.size();

	// Now, find the path of the optimal polygon or line.

	FloatPointPair path_point_j;
	PolygonVertex current_vertex;
	optimal_segments.clear();
	int i, j, k, index_2;

	// Vertices are in the order: i, j, k.
	if (closed)
	{
		j = potential_segs[0].prev_vertex;
		index_2 = j - start;
		if (index_2 < 0)
		{
			index_2 += n;
		}
		i = potential_segs[index_2].prev_vertex;
		k = potential_segs[0].vertex;
	}
	else {  // Open path.
		if (0 != start)
		{
			throw std::runtime_error("Non-zero start for open path.\n");
			return false;
		}
		j = potential_segs[n - 1].vertex;
		i = potential_segs[n - 1].prev_vertex;
		k = potential_segs[n - 1].vertex;
		index_2 = j;
	}
	path_point_j.x = path_points[j].x;
	path_point_j.y = path_points[j].y;
#ifdef VERTEX_ADJUSTMENT
	if (false == potential_segs[index_2].meeting_point)
	{
		path_point_j = AdjustVertex(i, j, k);
	}
#endif

	current_vertex.point.x = path_point_j.x;
	current_vertex.point.y = path_point_j.y;
	current_vertex.meeting_point = potential_segs[index_2].meeting_point;
	optimal_segments.insert(optimal_segments.begin(), current_vertex);

	superpixel->AddVertex(superpixel->XY2Pos(path_points[j]));

	k = j;
	j = potential_segs[index_2].prev_vertex;
	index_2 = j - start;
	if (index_2 < 0)
	{
		index_2 += n;
	}
	i = potential_segs[index_2].prev_vertex;
	while (index_2 >= 0)
	{
		path_point_j.x = path_points[j].x;
		path_point_j.y = path_points[j].y;
#ifdef VERTEX_ADJUSTMENT
		if (false == potential_segs[index_2].meeting_point)
		{
			path_point_j = AdjustVertex(i, j, k);
		}
#endif
		if ((abs(path_point_j.x - optimal_segments[0].point.x) > EFFECTIVE_ZERO) || (abs(path_point_j.y - optimal_segments[0].point.y) > EFFECTIVE_ZERO))
		{
			current_vertex.point.x = path_point_j.x;
			current_vertex.point.y = path_point_j.y;
			current_vertex.meeting_point = potential_segs[index_2].meeting_point;
			optimal_segments.insert(optimal_segments.begin(), current_vertex);
		}
		else {
			optimal_segments[0].meeting_point = (optimal_segments[0].meeting_point || potential_segs[index_2].meeting_point);
		}
		superpixel->AddVertex(superpixel->XY2Pos(path_points[j]));
		if (0 == index_2)
		{
			break;
		}
		k = j;
		j = potential_segs[index_2].prev_vertex;
		index_2 = j - start;
		if (index_2 < 0)
		{
			index_2 += n;
		}
		i = potential_segs[index_2].prev_vertex;
	}
	n = optimal_segments.size();
	if ((abs(optimal_segments[n - 1].point.x - optimal_segments[0].point.x) < EFFECTIVE_ZERO) && (abs(optimal_segments[n - 1].point.y - optimal_segments[0].point.y) < EFFECTIVE_ZERO))
	{
		optimal_segments.pop_back();
	}
	if (false == closed) // Extend endpoints.
	{
		SPixelData* processed_pixdat = superpixel->GetWorkspace()->GetPixeldata();
		FloatPointPair local_p0;
		FloatPointPair local_p1;
		PolygonVertex new_point;
		std::vector<PolygonVertex>::iterator seg_it;
		if (begin_endpoint)
		{
			seg_it = optimal_segments.begin();
			local_p0 = seg_it->point;
			new_point.meeting_point = false;
			if (optimal_segments.size() > 1)
			{
				seg_it++;
				local_p1 = seg_it->point;
			}
			else {
				local_p1 = local_p0;
				local_p1.x += 1.0;
			}
			new_point.point = processed_pixdat->AxialExtent(local_p0, local_p1, superpixel->GetIdentifier());
			optimal_segments.insert(optimal_segments.begin(), new_point);
		}
		if (end_endpoint)
		{
			if (optimal_segments.size() > 1)
			{
				seg_it = optimal_segments.end();
				seg_it--;
				local_p0 = seg_it->point;
				seg_it--;
				local_p1 = seg_it->point;
				new_point.meeting_point = false;
				new_point.point = processed_pixdat->AxialExtent(local_p0, local_p1, superpixel->GetIdentifier());
				optimal_segments.push_back(new_point);
			}
		}
	}
	return true;
}

bool Path::CalcSegmentPenalties()
{
	int n = potential_segs.size();
	if (n < 2)
	{
		return false;
	}
	int start = potential_segs[0].vertex;  // Offset for first vertex.
	potential_segs[0].prev_vertex = start;
	potential_segs[0].num_segments = 0;
	potential_segs[0].cumulative_penalty = 0;

	for (int counter = 1; counter < n; counter++)
	{
		potential_segs[counter].prev_vertex = potential_segs[counter].vertex;
		potential_segs[counter].num_segments = -1;
		potential_segs[counter].cumulative_penalty = 0;
	}
	for (int index_1 = 0; index_1 < n; index_1++)  // The first point in the potential segment.
	{
		std::vector<Segment> s;
		s = potential_segs[index_1].s;
		for (int k = 0; k < s.size(); k++)
		{
			int index_2 = s[k].j - start;  // The second point in a potential segment.
			if (index_2 < 0)
			{
				index_2 += n;
			}

			if ((potential_segs[index_2].num_segments <= 0) || (potential_segs[index_2].num_segments > (potential_segs[index_1].num_segments + 1)))
			{
				potential_segs[index_2].num_segments = potential_segs[index_1].num_segments + 1;
				potential_segs[index_2].prev_vertex = potential_segs[index_1].vertex;
				potential_segs[index_2].cumulative_penalty = potential_segs[index_1].cumulative_penalty + potential_segs[index_1].s[k].penalty;
			}
			else if (potential_segs[index_2].num_segments == (potential_segs[index_1].num_segments + 1))
			{
				if ((potential_segs[index_1].cumulative_penalty + potential_segs[index_1].s[k].penalty) < potential_segs[index_2].cumulative_penalty)
				{
					potential_segs[index_2].num_segments = potential_segs[index_1].num_segments + 1;
					potential_segs[index_2].prev_vertex = potential_segs[index_1].vertex;
					potential_segs[index_2].cumulative_penalty = potential_segs[index_1].cumulative_penalty + potential_segs[index_1].s[k].penalty;
				}
			}
		}
	}
	return true;
}

bool Path::Curve(bool glitch3)
{
	int n = optimal_segments.size();
	optimal_curve.clear();
	Corner local_corner;
	PolygonVertex p0, p1, p2;
	bool extend_p0, extend_p2;

	if (closed)
	{
		for (int i = 0; i < n; ++i)
		{
			if (i > 0)
			{
				p0 = optimal_segments[i - 1];
			}
			else {
				p0 = optimal_segments[n - 1];
			}
			p1 = optimal_segments[i];
			if (i < (n - 1))
			{
				p2 = optimal_segments[i + 1];
			}
			else {
				p2 = optimal_segments[0];
			}
			if (p1.meeting_point)
			{
				local_corner.smooth = false;
				local_corner.p0.x = (p0.point.x + p1.point.x) / 2.0;
				local_corner.p0.y = (p0.point.y + p1.point.y) / 2.0;
				local_corner.p1.x = (p1.point.x + p2.point.x) / 2.0;
				local_corner.p1.y = (p1.point.y + p2.point.y) / 2.0;
				local_corner.c0.x = p1.point.x;
				local_corner.c0.y = p1.point.y;
			}
			else {
				local_corner = CalculateCorner(p0.point, p1.point, p2.point);
			}
			optimal_curve.push_back(local_corner);
		}
	}
	else {  // Not closed, so line segments used for painting.
		if (n > 2)
		{
			for (int i = 1; i < (n - 1); ++i)
			{
				if (1 == i)
				{
					extend_p0 = true;
				}
				else {
					extend_p0 = false;
				}
				if ((n - 2) == i)
				{
					extend_p2 = true;
				}
				else {
					extend_p2 = false;
				}
				p0 = optimal_segments[i - 1];
				p1 = optimal_segments[i];
				p2 = optimal_segments[i + 1];
				local_corner = CalculateCorner(p0.point, p1.point, p2.point, 1.35, extend_p0, extend_p2, glitch3);
				optimal_curve.push_back(local_corner);
			}
		}
		else if (n > 1)
		{
			local_corner.p0 = optimal_segments[0].point;
			local_corner.p1 = optimal_segments[1].point;
			local_corner.c0.x = (optimal_segments[1].point.x + optimal_segments[0].point.x) / 2.0;
			local_corner.c0.y = (optimal_segments[1].point.y + optimal_segments[0].point.y) / 2.0;
			local_corner.c1 = local_corner.c0;
			local_corner.smooth = false;
			optimal_curve.push_back(local_corner);
		}
	}
	return true;
}

bool Path::InsertNext(Path* n)
{
	Path* current_next = next;
	next = n;
	n->SetPrev(this);
	n->SetNext(current_next);
	return true;
}

bool Path::InsertPrev(Path* p)
{
	Path* current_prev = prev;
	prev = p;
	p->SetNext(this);
	p->SetPrev(current_prev);
	return true;
}

bool Path::InsertSubpath(Path* p)
{
	if (NULL != subpath)
	{
		p->SetPrev(subpath->GetTail());
		p->SetNext(NULL);
		subpath->GetTail()->SetNext(p);
	}
	else {
		subpath = p;
		p->SetNext(NULL);
		p->SetPrev(NULL);
	}
	return true;
}

bool Path::SetNext(Path* n)
{
	next = n;
	return true;
}

bool Path::SetPrev(Path* p)
{
	prev = p;
	return true;
}

bool Path::MoveSuperPixel(SuperPixel* newSP)
{
	Path* old_prev = prev;
	Path* old_next = next;

	if (this == superpixel->GetPathHead())
	{
		superpixel->SetPathHead(old_next);
	}

	if (NULL != old_prev)
	{
		old_prev->next = old_next;
	}
	if (NULL != old_next)
	{
		old_next->prev = old_prev;
	}

	next = NULL;
	if (NULL == newSP->GetPathHead())
	{
		newSP->SetPathHead(this);
		prev = NULL;
	}
	else {
		newSP->SetPathTail();
		prev = newSP->GetPathTail();
		prev->next = this;
	}

	superpixel->SetPathTail();
	superpixel = newSP;
	superpixel->SetPathTail();
	return true;
}

Corner Path::CalculateCorner(FloatPointPair p0, FloatPointPair p1, FloatPointPair p2, float alpha_max, bool extend_p0, bool extend_p2, bool glitch3)
{
	// Input are the three points making up the corner: p0, p1, p2.
	// Ordinarily, the corner is based on the midpoints between p0 - p1 and p1 - p2.  But extend_p0 adn extend_p2 makes the actual p0 or p2 be used.
	Corner ret;
	// ret.p0 is the first midpoint.
	// ret.p1 is the second midpoint.
	// ret.c0 is the first control point, unless this is a corner, in which case it is the corner point.
	// ret.c1 is the second control point, if there is one.
	// ret.smooth indicates whether this is a curve (using control points) or a corner.
	// glitch3 causes some painted shapes to start with a radius of zero when they really shouldn't.
	FloatPointPair Mid0, Mid1, temp_point, intersection_point, closest_corner;
	float min_dist;

	// Extended points (at the ends) are brought in by one pixel so that it is less likely that the
	// radius calculation for the point will prematurely run into an edge.
	if (extend_p0)
	{
		float dx = p0.x - p1.x;
		float dy = p0.y - p1.y;
		float dist = sqrt(dx * dx + dy * dy);
		if (dist > 5)
		{
			Mid0.x = p1.x * (1.0 / dist) + p0.x * ((dist - 1.0) / dist);
			Mid0.y = p1.y * (1.0 / dist) + p0.y * ((dist - 1.0) / dist);
		}
		else {
			Mid0 = p0;
		}
	}
	else {
		Mid0.x = (p0.x + p1.x) / 2.0;
		Mid0.y = (p0.y + p1.y) / 2.0;
	}
	if (extend_p2)
	{
		float dx = p2.x - p1.x;
		float dy = p2.y - p1.y;
		float dist = sqrt(dx * dx + dy * dy);
		if (dist > 5)
		{
			Mid1.x = p1.x * (1.0 / dist) + p2.x * ((dist - 1.0) / dist);
			Mid1.y = p1.y * (1.0 / dist) + p2.y * ((dist - 1.0) / dist);
		}
		else {
			Mid1 = p2;
		}
	}
	else {
		Mid1.x = (p1.x + p2.x) / 2.0;
		Mid1.y = (p1.y + p2.y) / 2.0;
	}

	ret.p0.x = Mid0.x; // Set midpoints to the return value.
	ret.p0.y = Mid0.y;
	ret.p1.x = Mid1.x;
	ret.p1.y = Mid1.y;
	ret.c0.x = p1.x;
	ret.c0.y = p1.y;
	if (glitch3)
	{
		//ret.c1.x = p2.x;
		//ret.c1.y = p2.y;
	}
	else {
		ret.c1.x = p1.x;
		ret.c1.y = p1.y;
	}


	float mid_dx, int_dx;
	float mid_dy, int_dy;
	float inv_d02;
	float dist;

	mid_dx = Mid1.x - Mid0.x;
	mid_dy = Mid1.y - Mid0.y;
	if (((abs(mid_dx) + abs(mid_dy)) <= EFFECTIVE_ZERO) || // Midpoints are at the same point.
		((abs(p1.x - p0.x) <= EFFECTIVE_ZERO) && (abs(p1.y - p0.y) <= EFFECTIVE_ZERO)) || // p0 and p1 are at the same point.
		((abs(p2.x - p1.x) <= EFFECTIVE_ZERO) && (abs(p2.y - p1.y) <= EFFECTIVE_ZERO))) // p1 and p2 are at the same point.
	{
		ret.smooth = false;
		return ret;
	}

	// Try various possible intersection points on the p0-p1 line.
	// Start with the actual midpoint.
	intersection_point.x = Mid0.x;
	intersection_point.y = Mid0.y;

	inv_d02 = 1.0 / sqrt(mid_dx * mid_dx + mid_dy * mid_dy);  // Should have positive denominator, given earlier test of mid_dx and mid_dy.
	min_dist = abs(mid_dx * (p1.y - intersection_point.y) - mid_dy * (p1.x - intersection_point.x)) * inv_d02;  // Distance from first midpoint to p1.
	// Need to determine whether this intersection_point intersects with the unit square centered on p1.
	// If any parallel line passing through a corner is farther from the center than this, then then line intersects the square.
	bool intersects = false;
	for (int i = 0; (i < 4) && (false == intersects); i++)
	{
		if ((i & 1) == 1)
		{
			temp_point.x = p1.x - SMOOTH_BOX_SIZE;
		}
		else {
			temp_point.x = p1.x + SMOOTH_BOX_SIZE;
		}
		if ((i & 2) == 2)
		{
			temp_point.y = p1.y - SMOOTH_BOX_SIZE;
		}
		else {
			temp_point.y = p1.y + SMOOTH_BOX_SIZE;
		}
		dist = abs(mid_dx * (p1.y - temp_point.y) - mid_dy * (p1.x - temp_point.x)) * inv_d02; // Distance from p1 to the point on the parallel line passing through the corner.
		if (dist > min_dist)
		{
			intersects = true;  // I think this means we can just bail, as alpha will equal zero.  For now, don't make this optimization.
		}
	}
	if (false == intersects)
	{
		// Now, if no intersection, try points around p1 that are corners of a 1 unit square to see which is closest to midpoint line.
		for (int i = 0; i < 4; i++)
		{
			if ((i & 1) == 1)
			{
				temp_point.x = p1.x - SMOOTH_BOX_SIZE;
			}
			else {
				temp_point.x = p1.x + SMOOTH_BOX_SIZE;
			}
			if ((i & 2) == 2)
			{
				temp_point.y = p1.y - SMOOTH_BOX_SIZE;
			}
			else {
				temp_point.y = p1.y + SMOOTH_BOX_SIZE;
			}
			dist = abs(mid_dx * (Mid0.y - temp_point.y) - mid_dy * (Mid0.x - temp_point.x)) * inv_d02;  // Distance to the line passing through the midpoint Mid0.
			if ((i == 0) || (dist < min_dist))
			{
				closest_corner.x = temp_point.x;
				closest_corner.y = temp_point.y;
				min_dist = dist;
			}
		}
		// Calculate intersection_point from line parallel to Mid0-Mid1 but passing through corner.
		if (abs(mid_dx) > EFFECTIVE_ZERO)
		{
			float mid_beta = mid_dy / mid_dx;
			float line_alpha = closest_corner.y - (mid_beta * closest_corner.x);
			if (abs(p1.x - p0.x) > EFFECTIVE_ZERO)
			{
				float p0p1_beta = (p1.y - p0.y) / (p1.x - p0.x);
				float p0p1_alpha = p1.y - p0p1_beta * p1.x;
				intersection_point.x = (line_alpha - p0p1_alpha) / (p0p1_beta - mid_beta);  // Should not have near-zero denominator, since first intersection test was negative.
				intersection_point.y = line_alpha + mid_beta * intersection_point.x;
			}
			else { // Vertical p0p1
				intersection_point.y = line_alpha + mid_beta * p1.x;
				intersection_point.x = p1.x;
			}
		}
		else { // Vertical midpoint line
			float p0p1_beta = (p1.y - p0.y) / (p1.x - p0.x);  // Should not have near-zero denominator, since midpoint line is vertical.
			float p0p1_alpha = p1.y - p0p1_beta * p1.x;
			intersection_point.x = closest_corner.x;
			intersection_point.y = p0p1_alpha + p0p1_beta * intersection_point.x;
		}

	}
	// Now, intersection_point is where the line parallel to the midpoint line intersects the Mid0-p1 line.

	float p1_mid_dx = (p1.x - Mid0.x); // dx, dy for line from Mid0 to p1.
	float p1_mid_dy = (p1.y - Mid0.y);
	int_dx = (intersection_point.x - Mid0.x); // int_dx, int_dy for line from Mid0 to intersection_point.
	int_dy = (intersection_point.y - Mid0.y);
	float alpha = (4.0 / 3.0) * (sqrt(int_dx * int_dx + int_dy * int_dy) / sqrt(p1_mid_dx * p1_mid_dx + p1_mid_dy * p1_mid_dy));
	if ((alpha > alpha_max) || (alpha < CORNER_ALPHA_MIN))
	{
		ret.smooth = false;
		return ret;
	}
	else {
		ret.smooth = true;
		ret.c0.x = intersection_point.x;
		ret.c0.y = intersection_point.y;

		// Calculate intersection_point from line parallel to Mid0-Mid1 but passing through first intersection point.
		if (abs(mid_dx) > EFFECTIVE_ZERO)
		{
			float mid_beta = mid_dy / mid_dx;
			float line_alpha = intersection_point.y - (mid_beta * intersection_point.x);
			if (abs(p2.x - p1.x) > EFFECTIVE_ZERO)
			{
				float p1p2_beta = (p2.y - p1.y) / (p2.x - p1.x);
				float p1p2_alpha = p2.y - p1p2_beta * p2.x;
				ret.c1.x = (line_alpha - p1p2_alpha) / (p1p2_beta - mid_beta);  // Should not have near-zero denominator, since first intersection test was negative.
				ret.c1.y = line_alpha + mid_beta * ret.c1.x;
			}
			else { // Vertical p1p2
				ret.c1.y = line_alpha + mid_beta * p1.x;
				ret.c1.x = p1.x;
			}
		}
		else { // Vertical midpoint line
			float p1p2_beta = (p2.y - p1.y) / (p2.x - p1.x);  // Should not have near-zero denominator, since midpoint line is vertical.
			float p1p2_alpha = p2.y - p1p2_beta * p2.x;
			ret.c1.x = intersection_point.x;
			ret.c1.y = p1p2_alpha + p1p2_beta * ret.c1.x;
		}
	}
	return ret;
}

