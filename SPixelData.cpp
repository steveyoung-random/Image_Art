// Copyright (c) 2023-2025 Steve Young
// Licensed under the MIT License

#include "SPixelData.h"
#include "SuperPixel.h"
#include <cstring>

SPixelData::SPixelData(int w, int h)
{
	width = w;
	height = h;
	data = (int*)calloc(static_cast<size_t>(w) * h, sizeof(int));
	if ((NULL == data))
	{
		throw (std::runtime_error("Unable to allocate memory for SPixelData data.\n"));
	}
#ifdef USE_CUDA
	c_device_data = IntArray(w, h, true, 0);
#endif
	Reset();
}
SPixelData::SPixelData(const SPixelData& t)
{
	width = t.width;
	height = t.height;
	data = (int*)calloc(static_cast<size_t>(width) * height, sizeof(int));
	if ((NULL == data))
	{
		throw (std::runtime_error("Unable to allocate memory for SPixelData data.\n"));
	}
	int max_pos = width * height;
	for (int pos = 0; pos < max_pos; pos++)
	{
		data[pos] = t.data[pos];
	}
#ifdef USE_CUDA
	c_device_data = IntArray(width, height, false);
	if (!CopyFromHost(data, width * height, c_device_data))
	{
		throw std::runtime_error("Failed to copy SPixelData to c_device_data.\n");
	}
#endif
	meeting_points = t.meeting_points;
}

SPixelData::~SPixelData()
{
	if (NULL != data)
	{
		free(data);
	}
#ifdef USE_CUDA
	FreeIntArray(c_device_data);
#endif
	meeting_points.clear();
}

int SPixelData::GetPixel(int x, int y)
{
	if ((x >= 0) && (x < width) && (y >= 0) && (y < height))
	{
		return data[y * width + x];
	}
	return 0;
}

bool SPixelData::GetMeetingPoint(PointPair point)
{
	if (meeting_points.count(point.y * width + point.x) > 0)
	{
		return true;
	}
	return false;
}

int SPixelData::GetNeighbors(int identifier, PointPair point)
{
	int ret = 0;
	int x = point.x;
	int y = point.y;
	// Neighbor values around point o:
	//
	// 128   1   2
	//  64   o   4
	//  32  16   8
	//
	// Return value is the sum of values for surrounding pixels that have the value of identifier.
	//
	if (identifier != data[y * width + x])
	{
		return -1;
	}
	if (y > 0)
	{
		if (x > 0)
		{
			if (identifier == data[(y - 1) * width + (x - 1)])
			{
				ret += 128;
			}
		}
		if (identifier == data[(y - 1) * width + x])
		{
			ret += 1;
		}
		if (x < (width - 1))
		{
			if (identifier == data[(y - 1) * width + (x + 1)])
			{
				ret += 2;
			}
		}
	}
	if (x > 0)
	{
		if (identifier == data[y * width + (x - 1)])
		{
			ret += 64;
		}
	}
	if (x < (width - 1))
	{
		if (identifier == data[y * width + (x + 1)])
		{
			ret += 4;
		}
	}
	if (y < (height - 1))
	{
		if (x > 0)
		{
			if (identifier == data[(y + 1) * width + (x - 1)])
			{
				ret += 32;
			}
		}
		if (identifier == data[(y + 1) * width + x])
		{
			ret += 16;
		}
		if (x < (width - 1))
		{
			if (identifier == data[(y + 1) * width + (x + 1)])
			{
				ret += 8;
			}
		}
	}
	return ret;
}

bool SPixelData::SetMeetingPoint(PointPair point)
{
	meeting_points.insert(point.y * width + point.x);
	return true;
}
bool SPixelData::SetPixel(int x, int y, int value)
{
	if ((x >= 0) && (x < width) && (y >= 0) && (y < height))
	{
		data[y * width + x] = value;
		return true;
	}
	else {
//		throw (std::runtime_error("Attempt to set SPixelData value outside the boundaries.\n"));
		return false;
	}
}

bool SPixelData::Reset()
{
	int max_pos = width * height;
	memset(data, 0, max_pos * sizeof(int));
#ifdef USE_CUDA
	if (false == ResetIntArray(c_device_data, width, height, 0))
	{
		throw std::runtime_error("Failed to reset c_device_data in Reset.\n");
	}
	//c_device_data = IntArray(width, height, false);
	if (!CopyFromHost(data, max_pos, c_device_data))
	{
		throw std::runtime_error("Failed to copy SPixelData to c_device_data in Reset.\n");
	}
#endif
	meeting_points.clear();
	return true;
}

bool SPixelData::CopyData(SPixelData* sp)
{
	bool ret = true;
	if ((sp->height != height) || (sp->width != width))
	{
		throw std::runtime_error("Size mismatch for SPixelData Swap.\n");
		ret = false;
	}
	else {
		int max_pos = width * height;
		for (int pos = 0; pos < max_pos; pos++)
		{
			data[pos] = sp->data[pos];
		}
#ifdef USE_CUDA
		c_device_data = IntArray(width, height, false);
		if (!CopyFromHost(data, max_pos, c_device_data))
		{
			throw std::runtime_error("Failed to copy SPixelData to c_device_data in CopyData.\n");
			ret = false;
		}
#endif
	}
	return ret;
}

bool SPixelData::erode(int mode, int struct_size)
{
	// Mode: 0=square, 1=disc
	// Struct_size: Width/Height, or Diameter
	bool ret = false;
	SPixelData* local_data = new SPixelData(width, height);

	if (NULL == local_data)
	{
		throw (std::runtime_error("Unable to allocate memory for eroded SPixelData.\n"));
		return false;
	}
	if (struct_size < 3)
	{
		delete local_data;
		return false;
	}
	local_data->Reset();

	int windowsize = (struct_size - 1) / 2;
	bool set = true;
	bool last_set = true;
	int last_wx = -1 - windowsize;
	int* disc_chord;
	int identifier;

	// Calculate erosion
	if (0 == mode)		//Square
	{
		for (int j = 0; j < height; j++)
		{
			last_set = true;
			last_wx = -1 - windowsize;
			for (int i = 0; i < width; i++)
			{
				identifier = GetPixel(i, j);
				if ((last_wx < (i - windowsize)) || (true == last_set)) // Not using earlier results.
				{
					set = true;
					for (int wy = (j - windowsize); wy <= (j + windowsize); wy++)
					{
						if ((wy >= 0) && (wy < height))
						{
							for (int wx = (i - windowsize); wx <= (i + windowsize); wx++)
							{
								if ((wx >= 0) && (wx < width))
								{
									if (GetPixel(wx, wy) != identifier)
									{
										set = false;
										last_set = set;
										if (last_wx < wx)
										{
											last_wx = wx;
										}
									}
								}
							}
						}
					}
				}
				else {
					set = last_set;
					int start_wx = last_wx + 1;
					for (int wy = (j - windowsize); wy <= (j + windowsize); wy++)
					{
						if ((wy >= 0) && (wy < height))
						{
							for (int wx = start_wx; wx <= (i + windowsize); wx++)
							{
								if ((wx >= 0) && (wx < width))
								{
									if (GetPixel(wx, wy) != identifier)
									{
										set = false;
										last_set = set;
										if (last_wx < wx)
										{
											last_wx = wx;
										}
									}
								}
							}
						}
					}
				}
				if (set)
				{
					local_data->SetPixel(i, j, identifier);
				}
			}
		}
	}
	else
		// Disc
	{
		disc_chord = (int*)malloc(sizeof(int) * windowsize + sizeof(int));
		if (NULL == disc_chord)
		{
			throw (std::runtime_error("Unable to allocate memory for disc_chord for erode.\n"));
		}
		for (int i = 0; i <= windowsize; i++)
		{
			disc_chord[i] = sqrt(windowsize * windowsize - i * i);
		}
		for (int j = 0; j < height; j++)
		{
			for (int i = 0; i < width; i++)
			{
				set = true;
				identifier = GetPixel(i, j);
				for (int wy = (j - windowsize); wy <= (j + windowsize); wy++)
				{
					if ((wy >= 0) && (wy < height))
					{
						int disc_index = abs(wy - j);
						for (int wx = (i - disc_chord[disc_index]); wx <= (i + disc_chord[disc_index]); wx++)
						{
							if ((wx >= 0) && (wx < width))
							{
								if (GetPixel(wx, wy) != identifier)
								{
									set = false;
								}
							}
						}
					}
				}
				if (set)
				{
					local_data->SetPixel(i, j, identifier);
				}
			}
		}
		free(disc_chord);
	}
	ret = this->CopyData(local_data);
	delete local_data;
	return ret;
}

bool SPixelData::dilate(SuperPixel* sp, int mode, int struct_size)
{
	// Mode: 0=square, 1=disc
	// Struct_size: Width/Height, or Diameter
	bool ret = false;
	if (sp == NULL)
	{
		throw(std::runtime_error("Null SuperPixel provided to dilate function.\n"));
		return false;
	}
	SPixelData* local_SP = new SPixelData(width, height);
	if (NULL == local_SP)
	{
		throw (std::runtime_error("Unable to allocate memory for dilated SPixelData.\n"));
		return false;
	}
	if (struct_size < 3)
	{
		delete local_SP;
		return false;
	}
	local_SP->Reset();

	int windowsize = (struct_size - 1) / 2;
	bool set = false;
	bool last_set = false;
	int last_wx = -1 - windowsize;
	int* disc_chord = NULL;
	SuperPixel* current = NULL;
	current = sp->GetHead();

	if (1 == mode)
	{
		disc_chord = (int*)malloc(sizeof(int) * windowsize + sizeof(int));
		if (NULL == disc_chord)
		{
			throw (std::runtime_error("Unable to allocate memory for disc_chord for dilate.\n"));
		}
		for (int ii = 0; ii <= windowsize; ii++)
		{
			disc_chord[ii] = sqrt(windowsize * windowsize - ii * ii);
		}
	}
	while (NULL != current)
	{
		int identifier = current->GetIdentifier();
		RectQuad bbox = current->GetWindow();

		// Calculate dilation
		if (0 == mode)		//Square
		{
			for (int j = bbox.y0 - windowsize; j < (bbox.y1 + windowsize); j++)
			{
				if ((j >= 0) && (j < height))
				{
					last_set = false;
					last_wx = -1 - windowsize;
					for (int i = bbox.x0 - windowsize; i < (bbox.x1 + windowsize); i++)
					{
						if ((i >= 0) && (i < width))
						{
							if ((last_wx < (i - windowsize)) || (false == last_set)) // Not using earlier results.
							{
								set = false;
								for (int wy = j - windowsize; wy <= (j + windowsize); wy++)
								{
									if ((wy >= bbox.y0) && (wy <= bbox.y1))
									{
										for (int wx = i - windowsize; wx <= (i + windowsize); wx++)
										{
											if ((wx >= bbox.x0) && (wx <= bbox.x1))
											{
												if (GetPixel(wx, wy) == identifier)
												{
													set = true;
													last_set = set;
													last_wx = wx;
												}
											}
										}
									}
								}
							}
							else {
								set = last_set;
								int start_wx = last_wx + 1;
								for (int wy = j - windowsize; wy <= (j + windowsize); wy++)
								{
									if ((wy >= bbox.y0) && (wy <= bbox.y1))
									{
										for (int wx = start_wx; wx <= (i + windowsize); wx++)
										{
											if (wx <= bbox.x1)
											{
												if (GetPixel(wx, wy) == identifier)
												{
													set = true;
													last_set = set;
													last_wx = wx;
												}
											}
										}
									}
								}
							}
							if (set)
							{
								local_SP->SetPixel(i, j, identifier);
							}
						}
					}
				}
			}
		}
		else
			// Disc
		{
			for (int j = (bbox.y0 - windowsize); j <= (bbox.y1 + windowsize); j++)
			{
				if ((j >= 0) && (j < height))
				{
					for (int i = (bbox.x0 - windowsize); i <= (bbox.x1 + windowsize); i++)
					{
						if ((i >= 0) && (i < width))
						{
							set = false;
							for (int wy = (j - windowsize); wy <= (j + windowsize); wy++)
							{
								if ((wy >= bbox.y0) && (wy <= bbox.y1))
								{
									int disc_index = abs(wy - j);
									for (int wx = (i - disc_chord[disc_index]); wx <= (i + disc_chord[disc_index]); wx++)
									{
										if ((wx >= bbox.x0) && (wx <= bbox.x1))
										{
											if (GetPixel(wx, wy) == identifier)
											{
												set = true;
											}
										}
									}
								}
							}
							if (set)
							{
								local_SP->SetPixel(i, j, identifier);
							}
						}
					}
				}
			}
		}
		current = current->GetNext();
	}
	if (1 == mode)
	{
		free(disc_chord);
	}
	ret = this->CopyData(local_SP);
	delete local_SP;
	return ret;
}

bool SPixelData::dilate_erode(SuperPixel* sp, bool isdilate, int mode, int struct_size)
{

	// Mode: 0=square, 1=disc
	// Struct_size: Width/Height, or Diameter
	bool ret = false;
	if (sp == NULL)
	{
		throw(std::runtime_error("Null SuperPixel provided to dilate_erode function.\n"));
		return false;
	}
	SPixelData* local_SP = new SPixelData(width, height);
	if (NULL == local_SP)
	{
		throw (std::runtime_error("Unable to allocate memory for dilated or eroded SPixelData.\n"));
		return false;
	}
	if (struct_size < 3)
	{
		delete local_SP;
		return false;
	}
	local_SP->Reset();

	int windowsize = (struct_size - 1) / 2;
	if ((windowsize > (width - 1)) || (windowsize > (height - 1)))
	{
		throw std::runtime_error("Windowsize larger than image in dilate_erode.\n");
		delete local_SP;
		return false;
	}

	int* disc_chord = NULL;
	SuperPixel* current = NULL;
	current = sp->GetHead();

	if (1 == mode)
	{
		disc_chord = (int*)malloc(sizeof(int) * windowsize + sizeof(int));
		if (NULL == disc_chord)
		{
			throw (std::runtime_error("Unable to allocate memory for disc_chord for dilate_erode.\n"));
		}
		for (int ii = 0; ii <= windowsize; ii++)
		{
			disc_chord[ii] = sqrt(windowsize * windowsize - ii * ii);
		}
	}

	while (NULL != current)
	{
		int identifier = current->GetIdentifier();
		RectQuad bbox = current->GetWindow();
		int zero_sum = 0;
		int one_sum = 0;

		// Calculate dilation or erosion
		if (0 == mode)		//Square
		{
			int i = bbox.x0;
			int j = bbox.y0;
			int wx;
			int wy;
			// Calculate histogram for first pixel (bbox.x0,bbox.y0).
			for (wy = bbox.y0 - windowsize; wy <= (bbox.y0 + windowsize); wy++)
			{
				for (wx = bbox.x0 - windowsize; wx <= (bbox.x0 + windowsize); wx++)
				{
					if (GetPixel(wx, wy) == identifier)
					{
						one_sum++;
					}
					else {
						zero_sum++;
					}
				}
			}

			bool direction = true;  // True for moving to the right, false for moving to the left.

			while ((i <= bbox.x1) && (j <= bbox.y1))
			{
				// Now do a horizontal strip.
				if ((direction && (i > bbox.x0)) || ((false == direction) && (i < bbox.x1))) // If i == bbox.x0 or bbox.x1, then the histogram has already been filled out, no update needed.
				{
					for (wy = j - windowsize; wy <= (j + windowsize); wy++)
					{
						if (direction)
						{
							wx = i - (windowsize + 1); // Remove from histogram moving to the right.
						}
						else {
							wx = i + (windowsize + 1); // Remove from histogram moving to the left.
						}
						if (GetPixel(wx, wy) == identifier)
						{
							one_sum--;
						}
						else {
							zero_sum--;
						}
						if (direction)
						{
							wx = i + windowsize; // Add to histogram moving to the right.
						}
						else {
							wx = i - windowsize; // Add to histogram moving to the left.
						}
						if (GetPixel(wx, wy) == identifier)
						{
							one_sum++;
						}
						else {
							zero_sum++;
						}
					}
				}
				// Get the maximum or minimum from the histogram.
				if (isdilate)
				{
					if (one_sum > 0)
					{
						local_SP->SetPixel(i, j, identifier);
					}
				}
				else {
					if (zero_sum == 0)
					{
						local_SP->SetPixel(i, j, identifier);
					}
				}
				if (direction)
				{
					i++;
				}
				else {
					i--;
				}

				if (((i < bbox.x0) && (false == direction)) || ((i > bbox.x1) && direction)) // Need to drop down and reverse directions.
				{
					j++;
					if (j <= bbox.y1)
					{
						direction = !direction;
						if (direction)
						{
							i++;
						}
						else {
							i--;
						}
						// Calcluate changes to histogram from dropping down one row.
						for (wx = (i - windowsize); wx <= (i + windowsize); wx++)
						{
							wy = j - (windowsize + 1); // Remove from histogram
							if (GetPixel(wx, wy) == identifier)
							{
								one_sum--;
							}
							else {
								zero_sum--;
							}
							wy = j + windowsize; // Add to histogram
							if (GetPixel(wx, wy) == identifier)
							{
								one_sum++;
							}
							else {
								zero_sum++;
							}
						}
					}
				}
			}
		}
		else if (1 == mode)  // Disc
		{
			int i, j, wx, wy, disc_index;

			// Calculate histogram for first pixel (bbox.x0,bbox.y0).
			for (wy = bbox.y0 - windowsize; wy <= (bbox.y0 + windowsize); wy++)
			{
				disc_index = abs(wy - bbox.y0);
				for (wx = bbox.x0 - disc_chord[disc_index]; wx <= (bbox.x0 + disc_chord[disc_index]); wx++)
				{
					if (GetPixel(wx, wy) == identifier)
					{
						one_sum++;
					}
					else {
						zero_sum++;
					}
				}
			}

			i = bbox.x0;
			j = bbox.y0;
			bool direction = true;  // True for moving to the right, false for moving to the left.

			while ((i <= bbox.x1) && (j <= bbox.y1))
			{
				// Now do a horizontal strip.
				if ((direction && (i > bbox.x0)) || ((false == direction) && (i < bbox.x1))) // If i == bbox.x0 or bbox.x1, then the histogram may already have been filled out, no update needed.
				{
					for (wy = j - windowsize; wy <= (j + windowsize); wy++)
					{
						disc_index = abs(wy - j);
						if (direction)
						{
							wx = (i + disc_chord[disc_index]);
							if (GetPixel(wx, wy) == identifier)
							{
								one_sum++;
							}
							else {
								zero_sum++;
							}
							wx = (i - disc_chord[disc_index]) - 1;
							if (GetPixel(wx, wy) == identifier)
							{
								one_sum--;
							}
							else {
								zero_sum--;
							}
						}
						else {
							wx = (i + disc_chord[disc_index]) + 1;
							if (GetPixel(wx, wy) == identifier)
							{
								one_sum--;
							}
							else {
								zero_sum--;
							}
							wx = (i - disc_chord[disc_index]);
							if (GetPixel(wx, wy) == identifier)
							{
								one_sum++;
							}
							else {
								zero_sum++;
							}
						}
					}
				}
				// Get the maximum or minimum from the histogram.
				if (isdilate)
				{
					if (one_sum > 0)
					{
						local_SP->SetPixel(i, j, identifier);
					}
				}
				else {
					if (zero_sum == 0)
					{
						local_SP->SetPixel(i, j, identifier);
					}
				}
				if (direction)
				{
					i++;
				}
				else {
					i--;
				}

				if (((i < bbox.x0) && (false == direction)) || ((i > bbox.x1) && direction)) // Need to drop down and reverse directions.
				{
					j++;
					if (j <= bbox.y1)
					{
						direction = !direction;
						if (direction)
						{
							i++;
						}
						else {
							i--;
						}
						// Calcluate changes to histogram from dropping down one row.
						for (wx = (i - windowsize); wx <= (i + windowsize); wx++)
						{
							disc_index = abs(wx - i);
							wy = (j - disc_chord[disc_index]) - 1;
							if (GetPixel(wx, wy) == identifier)
							{
								one_sum--;
							}
							else {
								zero_sum--;
							}
							wy = (j + disc_chord[disc_index]);
							if (GetPixel(wx, wy) == identifier)
							{
								one_sum++;
							}
							else {
								zero_sum++;
							}
						}
					}
				}
			}
		}
		current = current->GetNext();
	}
	if (1 == mode)
	{
		free(disc_chord);
	}
	ret = this->CopyData(local_SP);
	delete local_SP;
	return ret;
}

ImageData* SPixelData::GenerateImage(SuperPixel* sp, Color background)
{
	ImageData* ret = NULL;
	SuperPixel* current = NULL;
	Color black = { 0, 0, 0 };
	Color red = { 255, 0, 0 };
	//black.channel[0] = 0;
	//black.channel[1] = 0;
	//black.channel[2] = 0;
	//red.channel[0] = 255;
	//red.channel[1] = 0;
	//red.channel[2] = 0;
	ret = new ImageData(NULL, width, height, 3, false);
	if (NULL == ret)
	{
		throw (std::runtime_error("Unable to allocate memory for generated image in SPixelData.\n"));
	}

	for (int j = 0; j < height; j++)
	{
		for (int i = 0; i < width; i++)
		{
			int sp_index = GetPixel(i, j);
			if (sp_index > 0)
			{
				current = sp->GetByIdentifier(sp_index);
				if (NULL == current)
				{
					throw std::runtime_error("Unable to find SuperPixel that matches value in pixeldata.\n");
				}
				Color c = current->GetAveColor();
				ret->SetPixel(i, j, c);
			}
			else if (sp_index == -1) // Skeleton endpoint
			{
				ret->SetPixel(i, j, black);
			}
			else if (sp_index == -2) // Skeleton intersection
			{
				ret->SetPixel(i, j, red);
			}
			else {
				ret->SetPixel(i, j, background);
			}
		}
	}

	return ret;
}

int SPixelData::GetWidth()
{
	return width;
}

int SPixelData::GetHeight()
{
	return height;
}

float SPixelData::CalculateRadius(std::vector<Corner>::iterator curve_begin, std::vector<Corner>::iterator curve_end, int mask_value)
{
	FloatPointPair p0, p1, c0, c1;
	std::vector<Corner>::iterator it;
	float last_radius = 0.0;
	float ret = 0.0;
	int count = 0;
	for (it = curve_begin; it != curve_end; ++it)
	{
		p0 = it->p0;
		p1 = it->p1;
		c0 = it->c0;
		c1 = it->c1;
		if (curve_begin == it)
		{
			it->radius_p0 = RadiusTransverse(p0, c0, mask_value);
			ret = it->radius_p0;
			count++;
		}
		else {
			it->radius_p0 = last_radius;
		}
		it->radius_p1 = RadiusTransverse(p1, c1, mask_value);
		ret += it->radius_p1;
		count++;
		last_radius = it->radius_p1;
		if (false == it->smooth) // Corner at c0.
		{
			// Change c1 to be on a line with c0 that is parallel to p0p1.
			c1.x = c0.x + p1.x - p0.x;
			c1.y = c0.y + p1.y - p0.y;
			it->radius_c0 = RadiusTransverse(c0, c1, mask_value);
			ret += it->radius_c0;
			count++;
		}
		else { // Proper bezier curve.  Assume that the apex is at t=0.5, though it may not be.  Hopefully it is close.
			FloatPointPair N = { 0, 0 }; // Point where apex is estimated to be.
			N.x = 0.125 * (p0.x + p1.x) + 0.375 * (c0.x + c1.x);
			N.y = 0.125 * (p0.y + p1.y) + 0.375 * (c0.y + c1.y);
			// Change c1 to be on a line with c0 that is parallel to p0p1.
			c1.x = N.x + p1.x - p0.x;
			c1.y = N.y + p1.y - p0.y;
			it->radius_c0 = RadiusTransverse(N, c1, mask_value);
			ret += it->radius_c0;
			count++;
		}
	}
	if (count > 0)
	{
		ret = ret / (float)count;
	}
	return ret;
}

float SPixelData::RadiusTransverse(FloatPointPair p, FloatPointPair c, int mask_value)
{
	// p - the point to measure from.
	// c - the control point defining the slope of the line.
	// mask_value - the value for the region in the data.

	float radius = 0.0;
	float dx = c.x - p.x;
	float dy = c.y - p.y;
	bool edge = false;
	float x_pos = p.x;
	float y_pos = p.y;
	float x_neg = p.x;
	float y_neg = p.y;
	float slope;  // Slope of the perpendicular line (so 90 degrees from the p-c line).
	int value;

	if ((abs(dx) < EFFECTIVE_ZERO) && (abs(dy) < EFFECTIVE_ZERO))
	{
		return 0.0;
	}

	if (abs(dy) > abs(dx)) // Perpendicular line is dominantly in the x direction.
	{
		slope = -dx / dy;
		while (false == edge)
		{
			x_pos += 1.0;
			y_pos += slope;
			x_neg -= 1.0;
			y_neg -= slope;
			if ((x_pos >= 0) && (x_pos < width) && (y_pos >= 0) && (y_pos < height))
			{
				value = GetPixel((int)x_pos, (int)y_pos);
				if (mask_value != value)
				{
					edge = true;
				}
			}
			else {
				edge = true;
			}
			if ((x_neg >= 0) && (x_neg < width) && (y_neg >= 0) && (y_neg < height))
			{
				value = GetPixel((int)x_neg, (int)y_neg);
				if (mask_value != value)
				{
					edge = true;
				}
			}
			else {
				edge = true;
			}
		}
		x_pos -= 1.0;
		y_pos -= slope;
		dx = x_pos - p.x;
		dy = y_pos - p.y;
		radius = sqrt(dx * dx + dy * dy);
	}
	else {                 // Perpendicular line is dominantly in the y direction.
		slope = -dy / dx;
		while (false == edge)
		{
			y_pos += 1.0;
			x_pos += slope;
			y_neg -= 1.0;
			x_neg -= slope;
			if ((x_pos >= 0) && (x_pos < width) && (y_pos >= 0) && (y_pos < height))
			{
				value = GetPixel((int)x_pos, (int)y_pos);
				if (mask_value != value)
				{
					edge = true;
				}
			}
			else {
				edge = true;
			}
			if ((x_neg >= 0) && (x_neg < width) && (y_neg >= 0) && (y_neg < height))
			{
				value = GetPixel((int)x_neg, (int)y_neg);
				if (mask_value != value)
				{
					edge = true;
				}
			}
			else {
				edge = true;
			}
		}
		y_pos -= 1.0;
		x_pos -= slope;
		dx = x_pos - p.x;
		dy = y_pos - p.y;
		radius = sqrt(dx * dx + dy * dy);
	}
	return radius;
}

FloatPointPair SPixelData::AxialExtent(FloatPointPair p, FloatPointPair c, int mask_value)
{
	// p - the point to measure from.
	// c - the control point defining the slope of the line, in the opposite direction from where we will measure the radius.
	// mask_value - the value for the region in the data.

	FloatPointPair ret = p;
	float dx = c.x - p.x;
	float dy = c.y - p.y;
	bool edge = false;
	float x_pos = p.x;
	float y_pos = p.y;
	float x_neg = p.x;
	float y_neg = p.y;
	float slope;  // Slope of the axial line (same as the p-c line).
	int value;
	bool move_pos;

	if ((abs(dx) < EFFECTIVE_ZERO) && (abs(dy) < EFFECTIVE_ZERO))
	{
		return ret;
	}

	if (abs(dy) > abs(dx)) // slope is dominantly in the y direction.
	{
		slope = dx / dy;
		if (dy > 0.0)
		{
			move_pos = false;
		}
		else {
			move_pos = true;
		}

		while (false == edge)
		{
			if (move_pos)
			{
				x_pos += slope;
				y_pos += 1.0;
				if ((x_pos >= 0) && (x_pos < width) && (y_pos >= 0) && (y_pos < height))
				{
					value = GetPixel((int)x_pos, (int)y_pos);
					if (mask_value != value)
					{
						edge = true;
					}
				}
				else {
					edge = true;
				}
			}
			else {
				x_neg -= slope;
				y_neg -= 1.0;
				if ((x_neg >= 0) && (x_neg < width) && (y_neg >= 0) && (y_neg < height))
				{
					value = GetPixel((int)x_neg, (int)y_neg);
					if (mask_value != value)
					{
						edge = true;
					}
				}
				else {
					edge = true;
				}
			}
		}
		if (move_pos)
		{
			x_pos -= slope;
			y_pos -= 1.0;
			ret.x = x_pos;
			ret.y = y_pos;
		}
		else {
			x_neg += slope;
			y_neg += 1.0;
			ret.x = x_neg;
			ret.y = y_neg;
		}
	}
	else {                 // slope  is dominantly in the x direction.
		slope = dy / dx;
		if (dx > 0.0)
		{
			move_pos = false;
		}
		else {
			move_pos = true;
		}
		while (false == edge)
		{
			if (move_pos)
			{
				y_pos += slope;
				x_pos += 1.0;
				if ((x_pos >= 0) && (x_pos < width) && (y_pos >= 0) && (y_pos < height))
				{
					value = GetPixel((int)x_pos, (int)y_pos);
					if (mask_value != value)
					{
						edge = true;
					}
				}
				else {
					edge = true;
				}
			}
			else {
				y_neg -= slope;
				x_neg -= 1.0;
				if ((x_neg >= 0) && (x_neg < width) && (y_neg >= 0) && (y_neg < height))
				{
					value = GetPixel((int)x_neg, (int)y_neg);
					if (mask_value != value)
					{
						edge = true;
					}
				}
				else {
					edge = true;
				}
			}
		}
		if (move_pos)
		{
			y_pos -= slope;
			x_pos -= 1.0;
			ret.x = x_pos;
			ret.y = y_pos;
		}
		else {
			y_neg += slope;
			x_neg += 1.0;
			ret.x = x_neg;
			ret.y = y_neg;
		}
	}
	return ret;
}

bool SPixelData::FloodReplace(int p, int orig, int updated)
{
	// Flood fill replace.  Starting at point p, replace all contiguous points that are value orig with the new value updated.
	// diagonals indicates whether to flood fill on diagonals as well as adjacent pixels.

	int pos = p;
	if (orig != data[pos])
	{
		throw std::runtime_error("FloodReplace called on point not in original set.\n");
		return false;
	}
	std::vector<int> remaining_set;
	remaining_set.clear();
	remaining_set.push_back(pos);

	while (remaining_set.size() > 0)
	{
		pos = remaining_set.back();
		remaining_set.pop_back();
		if (orig == data[pos])
		{
			data[pos] = updated;
			int y = pos / width;
			int x = pos - (y * width);
			if (y > 0)
			{
				remaining_set.push_back(pos - width);
			}
			if (y < (height - 1))
			{
				remaining_set.push_back(pos + width);
			}
			if (x > 0)
			{
				remaining_set.push_back(pos - 1);
			}
			if (x < (width - 1))
			{
				remaining_set.push_back(pos + 1);
			}
		}
	}
	return true;
}

#ifdef USE_CUDA
bool SPixelData::SyncToDevice()
{
	bool ret = true;
	if (!CopyFromHost(data, width*height, c_device_data))
	{
		throw std::runtime_error("Failed to copy SPixelData to c_device_data in SyncToDevice.\n");
		ret = false;
	}
	return ret;
}
int* SPixelData::GetDeviceData()
{
	return c_device_data;
}
#endif

