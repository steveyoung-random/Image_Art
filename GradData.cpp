// Copyright (c) 2023-2024 Steve Young
// Licensed under the MIT License

#include "GradData.h"
#include "ImageData.h"

GradData::GradData(unsigned char* gradient, int w, int h)
{
	width = w;
	height = h;
	data = (unsigned char*)malloc(sizeof(unsigned char) * width * height);
	if (NULL == data)
	{
		throw (std::runtime_error("Unable to allocate memory for GradData.\n"));
	}
	int total_pix = width * height;
	for (int pos = 0; pos < total_pix; pos++)
	{
		data[pos] = gradient[pos];
	}
}

GradData::GradData(const GradData& t)
{
	width = t.width;
	height = t.height;
	data = (unsigned char*)malloc(sizeof(unsigned char) * width * height);
	if (NULL == data)
	{
		throw (std::runtime_error("Unable to allocate memory for GradData.\n"));
	}
	long total_pix = width * height;
	for (long pos = 0; pos < total_pix; pos++)
	{
		data[pos] = t.data[pos];
	}
}

GradData::GradData(std::string filename)
{
	int colorchannels;
	unsigned char* data_temp;
	data_temp = stbi_load(filename.c_str(), &width, &height, &colorchannels, 0);
	if (NULL == data_temp)
	{
		throw (std::runtime_error("Unable to load grayscale image.\n"));
	}
	data = (unsigned char*)malloc(sizeof(unsigned char) * width * height);
	if (NULL == data)
	{
		throw (std::runtime_error("Unable to allocate memory for GradData.\n"));
	}
	int total_pix = width * height;
	for (int pos = 0; pos < total_pix; pos++)
	{
		data[pos] = data_temp[pos*colorchannels];
	}
	stbi_image_free(data_temp);
}

GradData::~GradData()
{
	if (NULL != data)
	{
		free(data);
	}
}

unsigned char GradData::GetPixel(int x, int y)
{
	return data[x + y * width];
}

int GradData::GetWidth()
{
	return width;
}

int GradData::GetHeight()
{
	return height;
}

bool GradData::write_file(std::string filename)
{
	std::stringstream ss;
	ss << filename << ".png";
	std::string s = ss.str();

	if (0 == stbi_write_png(s.c_str(), width, height, 1, data, width))
	{
		throw std::runtime_error("Unable to write out grayscale image.\n");
	}
	return true;
}

GradData* GradData::gen_dilate_erode(bool isdilate, int mode, int struct_size)
{
	// Mode: 0=square, 1=disc
	// Struct_size: Width/Height, or Diameter

	unsigned char* output = NULL;
	GradData* ret = NULL;
	int x = GetWidth();
	int y = GetHeight();
	int hist[256];

	output = (unsigned char*)malloc(sizeof(unsigned char) * x * y);
	if (NULL == output)
	{
		throw (std::runtime_error("Unable to allocate memory for image in gen_dilate_erode.\n"));
		return NULL;
	}

	for (int i = 0; i < 256; i++)
	{
		hist[i] = 0;
	}

	int windowsize = (struct_size - 1) / 2;
	if ((windowsize > (x - 1)) || (windowsize > (y - 1)))
	{
		throw std::runtime_error("Windowsize larger than image in gen_dilate_erode.\n");
		return NULL;
	}
	int* disc_chord;

	if (0 == mode)		//Square
	{
		int i = 0;
		int j = 0;
		int wx;
		int wy;
		// Calculate histogram for first pixel (0,0).
		for (wy = 0; wy <= windowsize; wy++)
		{
			for (wx = 0; wx <= windowsize; wx++)
			{
				hist[GetPixel(wx, wy)]++;
			}
		}

		int min_or_max;
		long pos;
		bool direction = true;  // True for moving to the right, false for moving to the left.

		while ((i < x) && (j < y))
		{
			// Now do a horizontal strip.
			if ((direction && (i > 0)) || ((false == direction) && (i < (x - 1)))) // If i == 0 or x-1, then the histogram has already been filled out, no update needed.
			{
				for (wy = j - windowsize; wy <= (j + windowsize); wy++)
				{
					if ((wy >= 0) && (wy < y))
					{
						if (direction)
						{
							wx = i - (windowsize + 1); // Remove from histogram moving to the right.
						}
						else {
							wx = i + (windowsize + 1); // Remove from histogram moving to the left.
						}
						if ((wx >= 0) && (wx < x))
						{
							hist[GetPixel(wx, wy)]--;
						}
						if (direction)
						{
							wx = i + windowsize; // Add to histogram moving to the right.
						}
						else {
							wx = i - windowsize; // Add to histogram moving to the left.
						}
						if ((wx >= 0) && (wx < x))
						{
							hist[GetPixel(wx, wy)]++;
						}
					}
				}
			}
			// Get the maximum or minimum from the histogram.
			if (isdilate)
			{
				min_or_max = 255;
				while (0 == hist[min_or_max])
				{
					min_or_max--;
					if (min_or_max < 0)
					{
						throw std::runtime_error("Histogram empty in gen_dilate_erode.\n");
						return NULL;
					}
				}
			}
			else {
				min_or_max = 0;
				while (0 == hist[min_or_max])
				{
					min_or_max++;
					if (min_or_max > 255)
					{
						throw std::runtime_error("Histogram empty in gen_dilate_erode.\n");
						return NULL;
					}
				}
			}
			pos = j * x + i;
			output[pos] = min_or_max;
			if (direction)
			{
				i++;
			}
			else {
				i--;
			}

			if (((i < 0) && (false == direction)) || ((i == x) && direction)) // Need to drop down and reverse directions.
			{
				j++;
				if (j < y)
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
						if ((wx >= 0) && (wx < x))
						{
							wy = j - (windowsize + 1); // Remove from histogram
							if (wy >= 0)
							{
								hist[GetPixel(wx, wy)]--;
							}
							wy = j + windowsize; // Add to histogram
							if (wy < y)
							{
								hist[GetPixel(wx, wy)]++;
							}
						}
					}
				}
			}
		}
	}
	else // Disc
	{
		int i, j, wx, wy, disc_index;

		disc_chord = (int*)malloc(sizeof(int) * (windowsize + 1));
		if (NULL == disc_chord)
		{
			throw (std::runtime_error("Unable to allocate memory for disc_chord for dilate.\n"));
		}
		for (i = 0; i <= windowsize; i++)
		{
			disc_chord[i] = sqrt(windowsize * windowsize - i * i);
		}

		// Calculate histogram for first pixel (0,0).
		for (wy = 0; wy <= windowsize; wy++)
		{
			disc_index = abs(wy);
			for (wx = 0; wx <= disc_chord[disc_index]; wx++)
			{
				hist[GetPixel(wx, wy)]++;
			}
		}
		i = 0;
		j = 0;
		int min_or_max;
		long pos;
		bool direction = true;  // True for moving to the right, false for moving to the left.

		while ((i < x) && (j < y))
		{
			// Now do a horizontal strip.
			if ((direction && (i > 0)) || ((false == direction) && (i < (x - 1)))) // If i == 0 or x-1, then the histogram has already been filled out, no update needed.
			{
				for (wy = j - windowsize; wy <= (j + windowsize); wy++)
				{
					if ((wy >= 0) && (wy < y))
					{
						disc_index = abs(wy - j);
						if (direction)
						{
							wx = (i + disc_chord[disc_index]);
							if (wx < x)
							{
								hist[GetPixel(wx, wy)]++;
							}
							wx = (i - disc_chord[disc_index]) - 1;
							if (wx >= 0)
							{
								hist[GetPixel(wx, wy)]--;
							}
						}
						else {
							wx = (i + disc_chord[disc_index]) + 1;
							if (wx < x)
							{
								hist[GetPixel(wx, wy)]--;
							}
							wx = (i - disc_chord[disc_index]);
							if (wx >= 0)
							{
								hist[GetPixel(wx, wy)]++;
							}
						}
					}
				}
			}
			// Get the maximum or minimum from the histogram.
			if (isdilate)
			{
				min_or_max = 255;
				while (0 == hist[min_or_max])
				{
					min_or_max--;
					if (min_or_max < 0)
					{
						throw std::runtime_error("Histogram empty in gen_dilate_erode.\n");
						return NULL;
					}
				}
			}
			else {
				min_or_max = 0;
				while (0 == hist[min_or_max])
				{
					min_or_max++;
					if (min_or_max > 255)
					{
						throw std::runtime_error("Histogram empty in gen_dilate_erode.\n");
						return NULL;
					}
				}
			}
			pos = j * x + i;
			output[pos] = min_or_max;
			if (direction)
			{
				i++;
			}
			else {
				i--;
			}

			if (((i < 0) && (false == direction)) || ((i == x) && direction)) // Need to drop down and reverse directions.
			{
				j++;
				if (j < y)
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
						if ((wx >= 0) && (wx < x))
						{
							disc_index = abs(wx - i);
							wy = (j - disc_chord[disc_index]) - 1;
							if (wy >= 0)
							{
								hist[GetPixel(wx, wy)]--;
							}
							wy = (j + disc_chord[disc_index]);
							if (wy < y)
							{
								hist[GetPixel(wx, wy)]++;
							}
						}
					}
				}
			}
		}
		free(disc_chord);
	}
	ret = new GradData(output, x, y);
	if (NULL == ret)
	{
		throw std::runtime_error("Failed to create GradData object in gen_dilate2.\n");
	}
	free(output);
	return ret;
}

GradData* GradData::gen_edge(GradData* dilate, GradData* erode, int xdiv, int ydiv, float konst)
{
	unsigned char* edge;
	GradData* ret = NULL;
	int x = dilate->GetWidth();
	int y = dilate->GetHeight();
	edge = (unsigned char*)malloc(sizeof(unsigned char) * x * y);

	if ((NULL == edge))
	{
		throw (std::runtime_error("Unable to allocate memory for edge image.\n"));
		return NULL;
	}

	// Calculate edge
	for (int j = 0; j < y; j++)
	{
		for (int i = 0; i < x; i++)
		{
			long pos = j * x + i;
			edge[pos] = dilate->GetPixel(i, j) - erode->GetPixel(i, j);
		}
	}

	// Calculate edge

	if (konst > 0)
	{
		for (int j = 0; j < y; j++)
		{
			int l = ydiv / 2 + (j / ydiv) * ydiv;
			for (int i = 0; i < x; i++)
			{
				int k = xdiv / 2 + (i / xdiv) * xdiv;
				long pos = j * x + i;
				float delta = konst * sqrt((float)(i - k) * (i - k) / (float)(xdiv * xdiv) + (float)(j - l) * (j - l) / (float)(ydiv * ydiv));
				if ((delta + (float)edge[pos]) >= 255)
				{
					edge[pos] = 255;
				}
				else {
					edge[pos] = edge[pos] + delta;
				}

			}
		}
	}
	ret = new GradData(edge, x, y);
	if (NULL == ret)
	{
		throw std::runtime_error("Failed to create GradData object in gen_edge.\n");
	}
	free(edge);
	return ret;
}

ImageData* GradData::Gradient2Image(int mode)
{
	//  Mode:  0=Gradient (edge)
	//         1+ = Gray

		ImageData* ret = NULL;
		unsigned char* data = NULL;
		unsigned char value = 0;
		if (mode < 0)
		{
			return ret;
		}
		data = (unsigned char*)malloc(sizeof(unsigned char) * width * height * 3);
		if (NULL == data)
		{
			throw std::runtime_error("Failed to allocate memory for Gradient2Image.\n");
		}
		for (int j = 0; j < height; j++)
		{
			for (int i = 0; i < width; i++)
			{
				int pos = 3 * (i + j * width);
				if (0 == mode)
				{
					value = GetPixel(i, j);
					if (value > 51)
					{
						value = 0;
					}
					else {
						value = 255 - 5 * value;
					}
				}
				else 
				{
					value = GetPixel(i, j);
				}
				data[pos] = value;
				data[pos + 1] = value;
				data[pos + 2] = value;  // Linear burn: add both values, subtract 255 (floor at zero).
			}
		}
		ret = new ImageData(data, width, height, 3, false);
		if (NULL == ret)
		{
			throw std::runtime_error("Failed to create ImageData in Gradient2Image.\n");
		}
		free(data);
		return ret;
}

GradData* GradData::Preprocess_Gray(int num_steps, unsigned char steps, unsigned char modes, int structsize)
{
	// num_steps - how many steps in the processing (0 to 7)
	// structsize - the size of the structuring elements used in processing.
	// steps - each bit (starting at 0 and going to 7) is 1 for erosion, 0 for dilation
	// modes - each bit (starting at 0 and going to 7) is 1 for disc, 0 for square
	// 

	GradData* processed_gray = NULL;
	GradData* temp_gray = NULL;

	processed_gray = new GradData(*this);

	if (0 == num_steps)
	{
		return processed_gray;
	}

	for (int i = 0; i < num_steps; i++)
	{
		int mode = 0;
		if (modes & (1 << i))
		{
			mode = 1;
		}
		else {
			mode = 0;
		}
		if (steps & (1 << i))  // True if erode
		{
			//temp_gray = gen_erode(processed_gray, mode, structsize);
			temp_gray = processed_gray->gen_dilate_erode(false, mode, structsize);
			delete(processed_gray);
			processed_gray = temp_gray;
			temp_gray = NULL;
		}
		else {  // Otherwise, dilate
			//temp_gray = gen_dilate(processed_gray, mode, structsize);
			temp_gray = processed_gray->gen_dilate_erode(true, mode, structsize);
			delete(processed_gray);
			processed_gray = temp_gray;
			temp_gray = NULL;
		}
	}
	return processed_gray;
}

GradData* GradData::Generate_Gradient(int mode, int struct_size, int xdiv, int ydiv, float konst)
{
	GradData* dilate = NULL;
	GradData* erode = NULL;
	GradData* edge = NULL;


	dilate = gen_dilate_erode(true, mode, struct_size);
	erode = gen_dilate_erode(false, mode, struct_size);

	edge = gen_edge(dilate, erode, xdiv, ydiv, konst);

	delete dilate;
	delete erode;

	return edge;
}

