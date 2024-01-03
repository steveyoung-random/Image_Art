#include "Workspace.h"

// Copyright (c) 2023-2024 Steve Young
// Licensed under the MIT License

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

GradData* WorkSpace::gen_diff(ImageData* image1, ImageData* image2)
{
	return image1->gen_diff(image2);
}

GradData* WorkSpace::gen_gray(ImageData* image, int channel, int nchannel)
{
	return image->gen_gray(channel, nchannel);
}


GradData* WorkSpace::gen_dilate_erode(GradData* gray, bool isdilate, int mode, int struct_size)
{
	return gray->gen_dilate_erode(isdilate, mode, struct_size);
}

GradData* WorkSpace::gen_dilate(GradData* gray, int mode, int struct_size)
{
	// Mode: 0=square, 1=disc
	// Struct_size: Width/Height, or Diameter

	unsigned char* dilate = NULL;
	GradData* ret = NULL;
	int x = gray->GetWidth();
	int y = gray->GetHeight();

	dilate = (unsigned char*)malloc(sizeof(unsigned char) * x * y);
	if (NULL == dilate)
	{
		throw (std::runtime_error("Unable to allocate memory for dilated image.\n"));
		return NULL;
	}

	int windowsize = (struct_size - 1) / 2;
	int max = 0;
	int last_max = 0;
	int last_wx = -1 - windowsize;
	PointPair last_point;
	int* disc_chord;

	last_point.x = last_wx;
	last_point.y = last_wx;

	// Calculate dilation
	if (0 == mode)
		//Square
	{
		for (int j = 0; j < y; j++)
		{
			last_max = 0;
			last_wx = -1 - windowsize;
			for (int i = 0; i < x; i++)
			{
				if ((last_wx < (i - windowsize)) || (0 == last_max))
				{
					max = 0;
					for (int wy = j - windowsize; wy <= (j + windowsize); wy++)
					{
						if ((wy >= 0) && (wy < y))
						{
							for (int wx = i - windowsize; wx <= (i + windowsize); wx++)
							{
								if ((wx >= 0) && (wx < x))
								{
									if (gray->GetPixel(wx, wy) >= max)
									{
										max = gray->GetPixel(wx, wy);
										last_max = max;
										last_wx = wx;
									}
								}
							}
						}
					}
				}
				else {
					max = last_max;
					int start_wx = last_wx + 1;
					for (int wy = j - windowsize; wy <= (j + windowsize); wy++)
					{
						if ((wy >= 0) && (wy < y))
						{
							for (int wx = start_wx; wx <= (i + windowsize); wx++)
							{
								if (wx < x)
								{
									if (gray->GetPixel(wx, wy) >= max)
									{
										max = gray->GetPixel(wx, wy);
										last_max = max;
										last_wx = wx;
									}
								}
							}
						}
					}
				}
				long pos = j * x + i;
				dilate[pos] = max;
			}
		}
	}
	else
		// Disc
	{
		disc_chord = (int*)malloc(sizeof(int) * (windowsize + 1));
		if (NULL == disc_chord)
		{
			throw (std::runtime_error("Unable to allocate memory for disc_chord for dilate.\n"));
		}
		for (int i = 0; i <= windowsize; i++)
		{
			disc_chord[i] = sqrt(windowsize * windowsize - i * i);
		}
		for (int j = 0; j < y; j++)
		{
			last_max = 0;
			last_point.x = -1 - windowsize;
			last_point.y = last_point.x;
			for (int i = 0; i < x; i++)
			{
				// if (last_point.y < (j - windowsize));  // last point no longer in scope in y direction.
				// if (last_point.x < (i - disc_chord[abs(last_point.y - j)])); // last point not in scope in x direction.
				// if (0 == last_max); // There is no last point.
				// first wx is no further left than former right curve:
				//  old i = last_point.x - disc_chord[abs(last_point.y-j)]
				//  right curve = old i + disc_chord[disc_index]
				//  right curve = last_point.x - disc_chord[abs(last_point.y-j)] + disc_chord[disc_index]

				if ((last_point.y < (j - windowsize)) || (last_point.x < (i - disc_chord[abs(last_point.y - j)])) || (0 == last_max))
				{
					max = 0;
					for (int wy = j - windowsize; wy <= (j + windowsize); wy++)
					{
						if ((wy >= 0) && (wy < y))
						{
							int disc_index = abs(wy - j);
							for (int wx = i - disc_chord[disc_index]; wx <= (i + disc_chord[disc_index]); wx++)
							{
								if ((wx >= 0) && (wx < x))
								{
									if (gray->GetPixel(wx, wy) > max)
									{
										max = gray->GetPixel(wx, wy);
										last_max = max;
										last_point.x = wx;
										last_point.y = wy;
									}
								}
							}
						}
					}
				}
				else {
					max = last_max;
					int lx = last_point.x;
					int ly = last_point.y;
					for (int wy = j - windowsize; wy <= (j + windowsize); wy++)
					{
						if ((wy >= 0) && (wy < y))
						{
							int disc_index = abs(wy - j);
							int start_wx = lx - disc_chord[abs(ly - j)] + disc_chord[disc_index] + 1;
							start_wx = std::max(start_wx, (i - disc_chord[disc_index]));
							for (int wx = start_wx; wx <= (i + disc_chord[disc_index]); wx++)
							{
								if ((wx >= 0) && (wx < x))
								{
									if (gray->GetPixel(wx, wy) > max)
									{
										max = gray->GetPixel(wx, wy);
										last_max = max;
										last_point.x = wx;
										last_point.y = wy;
									}
								}
							}
						}
					}
				}
				long pos = j * x + i;
				dilate[pos] = max;
			}
		}
		free(disc_chord);
	}
	ret = new GradData(dilate, x, y);
	if (NULL == ret)
	{
		throw std::runtime_error("Failed to create GradData object in gen_dilate.\n");
	}
	free(dilate);
	return ret;
}

GradData* WorkSpace::gen_erode(GradData* gray, int mode, int struct_size)
{
	// Mode: 0=square, 1=disc
	// Struct_size: Width/Height, or Diameter
	unsigned char* erode = NULL;
	GradData* ret = NULL;
	int x = gray->GetWidth();
	int y = gray->GetHeight();
	erode = (unsigned char*)malloc(sizeof(unsigned char) * x * y);
	if (NULL == erode)
	{
		throw (std::runtime_error("Unable to allocate memory for eroded image.\n"));
		return NULL;
	}

	int windowsize = (struct_size - 1) / 2;
	int min = 255;
	int last_min = 255;
	int last_wx = -1 - windowsize;
	PointPair last_point;
	int* disc_chord;

	last_point.x = last_wx;
	last_point.y = last_wx;

	// Calculate erosion
	if (0 == mode)
		//Square
	{
		for (int j = 0; j < y; j++)
		{
			last_min = 255;
			last_wx = -1 - windowsize;
			for (int i = 0; i < x; i++)
			{
				if ((last_wx < (i - windowsize)) || (255 == last_min))
				{
					min = 255;
					for (int wy = j - windowsize; wy <= (j + windowsize); wy++)
					{
						if ((wy >= 0) && (wy < y))
						{
							for (int wx = i - windowsize; wx <= (i + windowsize); wx++)
							{
								if ((wx >= 0) && (wx < x))
								{
									if (gray->GetPixel(wx, wy) <= min)
									{
										min = gray->GetPixel(wx, wy);
										last_min = min;
										last_wx = wx;
									}
								}
							}
						}
					}
				}
				else {
					min = last_min;
					int start_wx = last_wx + 1;
					for (int wy = j - windowsize; wy <= (j + windowsize); wy++)
					{
						if ((wy >= 0) && (wy < y))
						{
							for (int wx = start_wx; wx <= (i + windowsize); wx++)
							{
								if (wx < x)
								{
									if (gray->GetPixel(wx, wy) <= min)
									{
										min = gray->GetPixel(wx, wy);
										last_min = min;
										last_wx = wx;
									}
								}
							}
						}
					}
				}
				long pos = j * x + i;
				erode[pos] = min;
			}
		}
	}
	else
		// Disc
	{
		disc_chord = (int*)malloc(sizeof(int) * (windowsize + 1));
		if (NULL == disc_chord)
		{
			throw (std::runtime_error("Unable to allocate memory for disc_chord for erode.\n"));
		}
		for (int i = 0; i <= windowsize; i++)
		{
			disc_chord[i] = sqrt(windowsize * windowsize - i * i);
		}
		for (int j = 0; j < y; j++)
		{
			last_min = 255;
			last_point.x = -1 - windowsize;
			last_point.y = last_point.x;
			for (int i = 0; i < x; i++)
			{
				if ((last_point.y < (j - windowsize)) || (last_point.x < (i - disc_chord[abs(last_point.y - j)])) || (255 == last_min))
				{
					min = 255;
					for (int wy = j - windowsize; wy <= (j + windowsize); wy++)
					{
						if ((wy >= 0) && (wy < y))
						{
							int disc_index = abs(wy - j);
							for (int wx = i - disc_chord[disc_index]; wx <= (i + disc_chord[disc_index]); wx++)
							{
								if ((wx >= 0) && (wx < x))
								{
									if (gray->GetPixel(wx, wy) < min)
									{
										min = gray->GetPixel(wx, wy);
										last_min = min;
										last_point.x = wx;
										last_point.y = wy;
									}
								}
							}
						}
					}
				}
				else {
					min = last_min;
					int lx = last_point.x;
					int ly = last_point.y;
					for (int wy = j - windowsize; wy <= (j + windowsize); wy++)
					{
						if ((wy >= 0) && (wy < y))
						{
							int disc_index = abs(wy - j);
							int start_wx = lx - disc_chord[abs(ly - j)] + disc_chord[disc_index] + 1;
							start_wx = std::max(start_wx, (i - disc_chord[disc_index]));
							for (int wx = start_wx; wx <= (i + disc_chord[disc_index]); wx++)
							{
								if ((wx >= 0) && (wx < x))
								{
									if (gray->GetPixel(wx, wy) < min)
									{
										min = gray->GetPixel(wx, wy);
										last_min = min;
										last_point.x = wx;
										last_point.y = wy;
									}
								}
							}
						}
					}
				}

				long pos = j * x + i;
				erode[pos] = min;
			}
		}
		free(disc_chord);
	}

	ret = new GradData(erode, x, y);
	if (NULL == ret)
	{
		throw std::runtime_error("Failed to create GradData object in gen_erode.\n");
	}
	free(erode);
	return ret;
}

GradData* WorkSpace::gen_edge(GradData* dilate, GradData* erode, int xdiv, int ydiv, float konst)
{
	return dilate->gen_edge(dilate, erode, xdiv, ydiv, konst);
}

int WorkSpace::GetXdiv()
{
	return xdiv;
}

int WorkSpace::GetYdiv()
{
	return ydiv;
}

SuperPixel* WorkSpace::GetHead()
{
	return list_head;
}

ImageData* WorkSpace::GetImage()
{
	return image;
}

GradData* WorkSpace::GetGray()
{
	return gray;
}

SPixelData* WorkSpace::GetPixeldata()
{
	if (NULL != processed_pixeldata)
	{
		return processed_pixeldata;
	}
	return pixeldata;
}

bool WorkSpace::SetXdiv(int xd)
{
	xdiv = xd;
	return true;
}

bool WorkSpace::SetYdiv(int yd)
{
	ydiv = yd;
	return true;
}



WorkSpace::WorkSpace(std::string filename, int channel, int nchannel, bool d)
{
	data = stbi_load(filename.c_str(), &width, &height, &colorchannels, 0);
	if (NULL != data)
	{
		std::cout << "X: " << width << " Y: " << height << " Colorchannels: " << colorchannels << "\n";
		bbox.x0 = 0;
		bbox.y0 = 0;
		bbox.x1 = width - 1;
		bbox.y1 = height - 1;
	}
	else {
		throw (std::runtime_error("Unable to load image.\n"));
	}

	image = new ImageData(data, width, height, colorchannels);
	if (NULL == image)
	{
		throw (std::runtime_error("Failed to create ImageData object.\n"));
	}
	gray = image->gen_gray(channel, nchannel);
	pixeldata = new SPixelData(width, height);
	diagonals = d;
	superpixel_set.clear();
}

WorkSpace::WorkSpace(std::string filename, std::string sp_filename, std::string gray_filename, std::string edge_filename, bool d)
{

	data = stbi_load(filename.c_str(), &width, &height, &colorchannels, 0);
	if (NULL != data)
	{
		std::cout << "X: " << width << " Y: " << height << " Colorchannels: " << colorchannels << "\n";
		bbox.x0 = 0;
		bbox.y0 = 0;
		bbox.x1 = width - 1;
		bbox.y1 = height - 1;
	}
	else {
		throw (std::runtime_error("Unable to load image.\n"));
	}

	image = new ImageData(data, width, height, colorchannels);
	if (NULL == image)
	{
		throw (std::runtime_error("Failed to create ImageData object.\n"));
	}

	if ("" == gray_filename)
	{
		gray = image->gen_gray();
	}
	else {
		gray = new GradData(gray_filename);
	}

	if ("" == edge_filename)
	{
		throw std::runtime_error("Failed to supply grayscale image filename.\n");
	}
	else {
		edge = new GradData(edge_filename);
	}

	pixeldata = new SPixelData(width, height);
	diagonals = d;
	list_head = NULL;
	list_tail = NULL;

	std::ifstream ifile;
	ifile.open(sp_filename, std::ios::binary);
	bool more = true;

	int identifier;
	PointPair seed;
	RectQuad box;
	Color color;
	int sp_width, size;

	while (more)
	{
		ifile.read((char*)(&identifier), sizeof(identifier));
		ifile.read((char*)(&seed.x), sizeof(seed.x));
		ifile.read((char*)(&seed.y), sizeof(seed.y));
		ifile.read((char*)(&box.x0), sizeof(box.x0));
		ifile.read((char*)(&box.y0), sizeof(box.y0));
		ifile.read((char*)(&box.x1), sizeof(box.x1));
		ifile.read((char*)(&box.y1), sizeof(box.y1));
		ifile.read((char*)(&color.channel[0]), sizeof(color.channel[0]));
		ifile.read((char*)(&color.channel[1]), sizeof(color.channel[1]));
		ifile.read((char*)(&color.channel[2]), sizeof(color.channel[2]));
		sp_width = box.x1 - box.x0 + 1;
		size = 0;
		for (int j = box.y0; j <= box.y1; j++)
		{
			int i = box.x0;
			int sequence;
			bool on = false;  // Each line starts with off.
			ifile.read((char*)(&sequence), sizeof(sequence));
			if (sequence > sp_width)
			{
				throw std::runtime_error("Error reading from SuperPixel file.  Too many leading zeros.\n");
				return;
			}
			i = i + sequence;  // No writing zeros, since other SuperPixels may write to those locations.
			while (i < (box.x1 + 1))
			{
				on = !on;
				ifile.read((char*)(&sequence), sizeof(sequence));
				if ((sequence + i - 1) > box.x1)
				{
					throw std::runtime_error("Error reading from SuperPixel file.  Too many items in sequence.\n");
					return;
				}
				if (on)
				{
					size = size + sequence;
					while (sequence > 0)
					{
						pixeldata->SetPixel(i, j, identifier);
						i++;
						sequence--;
					}
				}
				else {
					i = i + sequence;
				}
			}
			if (i > (box.x1 + 1))  // Numbers of sequence should add up to exactly length of line.
			{
				throw std::runtime_error("Error reading from SuperPixel file.  Items in sequence ran past end of line.\n");
				return;
			}
		}

		if (NULL == list_head)
		{
			list_head = new SuperPixel(identifier, edge, pixeldata, seed, NULL, NULL, this, SPType_Processed);
			InsertSPIntoIndex(list_head);
			list_tail = list_head;
		}
		else {
			list_tail = new SuperPixel(identifier, edge, pixeldata, seed, NULL, list_head->GetTail(), this, SPType_Processed);
			InsertSPIntoIndex(list_tail);
		}
		list_tail->SetWindow(box);
		list_tail->SetSize(size);
		list_tail->SetAveColor(color);
		// Test for end of file.
		more = (EOF != ifile.peek());
	}
	ifile.close();
}

WorkSpace::~WorkSpace()
{
	if (NULL != list_head)
	{
		current = list_head->GetHead();
		while (NULL != current->GetNext())
		{
			current = current->GetNext();
			delete current->GetPrev();
		}
		if (NULL != current)
		{
			delete current;
		}
	}
	if (NULL != skeleton_head)
	{
		current = skeleton_head->GetHead();
		while (NULL != current->GetNext())
		{
			current = current->GetNext();
			delete current->GetPrev();
		}
		if (NULL != current)
		{
			delete current;
		}
	}
	if (NULL != processed_head)
	{
		current = processed_head->GetHead();
		while (NULL != current->GetNext())
		{
			current = current->GetNext();
			delete current->GetPrev();
		}
		if (NULL != current)
		{
			delete current;
		}
	}
	if (gray)
	{
		delete(gray);
	}
	if (processed_gray)
	{
		delete(processed_gray);
	}
	if (dilate)
	{
		delete(dilate);
	}
	if (temp_gray)
	{
		delete(temp_gray);
	}
	if (erode)
	{
		delete(erode);
	}
	if (edge) {
		delete(edge);
	}
	if (data_revised)
	{
		delete(data_revised);
	}
	if (data_diff)
	{
		delete(data_diff);
	}
	if (pixeldata)
	{
		delete pixeldata;
	}
	if (skeleton_pixeldata)
	{
		delete skeleton_pixeldata;
	}
	if (processed_pixeldata)
	{
		delete processed_pixeldata;
	}
	if (image)
	{
		delete(image);
	}
	if (data)
	{
		stbi_image_free(data);
	}
	if (NULL != palette)
	{
		delete(palette);
	}
}

bool WorkSpace::Reset()
{
	return false;
}

bool WorkSpace::write_file(std::string filename)
{
	return false;
}

bool WorkSpace::Preprocess_Gray(int num_steps, unsigned char steps, unsigned char modes, int structsize)
{
	processed_gray = gray->Preprocess_Gray(num_steps, steps, modes, structsize);
	return true;
}

bool WorkSpace::Generate_Gradient(int mode, int struct_size, int xd, int yd, float konst)
{
	xdiv = xd;
	ydiv = yd;

	if (NULL != dilate)
	{
		delete(dilate);
	}
	dilate = processed_gray->gen_dilate_erode(true, mode, struct_size);

	if (NULL != erode)
	{
		delete(erode);
	}
	erode = processed_gray->gen_dilate_erode(false, mode, struct_size);

	if (NULL != edge)
	{
		delete(edge);
	}
	edge = gen_edge(dilate, erode, xdiv, ydiv, konst);
	return true;
}

bool WorkSpace::InitialSuperPixels(std::string seeds)
{
	std::ifstream ifile;

	if (NULL != list_head)
	{
		current = list_head->GetHead();
		while (NULL != current->GetNext())
		{
			current = current->GetNext();
			delete current->GetPrev();
		}
		if (NULL != current)
		{
			delete current;
		}
	}
	list_head = NULL;
	list_tail = NULL;
	int count = 0;

	if ("" != seeds)
	{
		PointPair seed;
		ifile.open(seeds, std::ios::in);
		if (ifile.is_open())
		{
			while (ifile >> seed.x >> seed.y)
			{
				count++;
				if (NULL == list_head)
				{
					list_head = new SuperPixel(count, edge, pixeldata, seed, NULL, NULL, this, SPType_Plain);
					InsertSPIntoIndex(list_head);
					list_tail = list_head;
				}
				else {
					list_tail = new SuperPixel(count, edge, pixeldata, seed, NULL, list_head->GetTail(), this, SPType_Plain);
					InsertSPIntoIndex(list_tail);
				}
			}
		}
		else {
			throw std::runtime_error("Failed to open seed file.\n");
			return false;
		}
		ifile.close();
	}
	else {
		for (int j = ydiv / 2; j < height; j = j + ydiv)
		{
			for (int i = xdiv / 2; i < width; i = i + xdiv)
			{
				count++;
				cell = new Cell(edge, pixeldata, i, j, xdiv, ydiv, 2);
				if (cell->FindSeed(NULL, 0, diagonals))
				{
					PointPair seed = cell->GetSeed();
					std::cout << "Seed: " << seed.x << ", " << seed.y << "\n";
					delete cell;
					if (NULL == list_head)
					{
						list_head = new SuperPixel(count, edge, pixeldata, seed, NULL, NULL, this, SPType_Plain);
						InsertSPIntoIndex(list_head);
						list_tail = list_head;
					}
					else {
						list_tail = new SuperPixel(count, edge, pixeldata, seed, NULL, list_head->GetTail(), this, SPType_Plain);
						InsertSPIntoIndex(list_tail);
					}
				}
				else {
					delete cell;
				}
			}
		}
	}
	return true;
}

bool WorkSpace::Watershed()
{
	pixeldata->Reset();
	current = list_head;
	while (NULL != current)
	{
		current->Reset();
		current = current->GetNext();
	}

	// Implement actual waterpixel algorithm.
	for (int level = 0; level < 256; level++)
	{
		bool done = false;
		std::cout << "Starting level: " << level << "\n";
		while (false == done)
		{
			done = true;
			current = list_head;
			while (NULL != current)
			{
				int grow_ret = current->Grow(level, false, true, bbox, NULL, 0, diagonals);
				if (grow_ret > 0)
				{
					done = false;
				}
				current = current->GetNext();
			}
		}
	}
	return true;
}

bool WorkSpace::SetAveColors()
{
	if (NULL != data_revised)
	{
		delete(data_revised);
	}
	data_revised = new ImageData(NULL, width, height, colorchannels);
	if (NULL == data_revised)
	{
		throw (std::runtime_error("Unable to allocate memory for revised image.\n"));
	}

	// Set seed color of each SuperPixel to average color.
	current = list_head;
	while (NULL != current)
	{
		Color ave;
		ave = current->SetAveColor(image);
		PointPair seed = current->GetSeed();

		std::cout << current->GetIdentifier() << ", " << seed.x << ", " << seed.y << ", " << current->GetSize() << ", " << current->GetAveError() << "\n";
		if (current->SetNeighbors())
		{
			std::set<int>* neighbors;
			neighbors = current->GetNeighbors();
			std::set<int>::iterator it;
			for (it = neighbors->begin(); it != neighbors->end(); ++it)
			{
				std::cout << *it << " ";
			}
			std::cout << "\n";
		}
		current = current->GetNext();
	}
	return true;
}

bool WorkSpace::CombineSuperPixels(float colormatch)
{
	bool done = false;
	while (false == done)
	{
		current = list_head;
		done = true;
		while (NULL != current)
		{
			Color ave1, ave2;
			SuperPixel* nsp;
			ave1 = current->GetAveColor();
			std::set<int>* neighbors;
			neighbors = current->GetNeighbors();
			for (std::set<int>::iterator iter = neighbors->begin(); iter != neighbors->end(); ++iter)
			{
				nsp = current->GetByIdentifier(*iter);
				if (NULL != nsp)
				{
					ave2 = nsp->GetAveColor();
					if (current->ColorDifference(ave1, ave2) < colormatch)
					{
						done = false;
						if (current->GetSize() >= nsp->GetSize())
						{
							current->Absorb(nsp, true, true);
							list_head = current->GetHead();
							list_tail = current->GetTail();
							break;
						}
					}
				}
			}
			current = current->GetNext();
		}
	}
	return true;
}

bool WorkSpace::Postprocess_SuperPixels(int num_steps, unsigned char steps, unsigned char modes, int structsize)
{
	// num_steps - how many steps in the processing (0 to 7)
	// windowsize - the size of the structuring elements used in processing.
	// steps - each bit (starting at 0 and going to 7) is 1 for erosion, 0 for dilation
	// modes - each bit (starting at 0 and going to 7) is 1 for disc, 0 for square
	// 

	if (0 == num_steps)
	{
		return true;
	}
	SuperPixel* processed_current = NULL;
	if (NULL == processed_pixeldata)
	{
		processed_pixeldata = new SPixelData(*pixeldata);
	}
	if (NULL == processed_head)
	{
		current = list_head;
		processed_current = new SuperPixel(*current, edge, processed_pixeldata, NULL, NULL, this, SPType_Processed);
		processed_head = processed_current;
		while (NULL != current->GetNext())
		{
			current = current->GetNext();
			processed_current = new SuperPixel(*current, edge, processed_pixeldata, NULL, processed_current, this, SPType_Processed);
		}
		processed_tail = processed_head->GetTail();
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
			if (false == processed_pixeldata->dilate_erode(processed_head, false, mode, structsize))
			{
				throw std::runtime_error("Error eroding SuperPixels.\n");
			}
		}
		else {  // Otherwise, dilate
			if (false == processed_pixeldata->dilate_erode(processed_head, true, mode, structsize))
			{
				throw std::runtime_error("Error dilating SuperPixels.\n");
			}
		}
	}
	return true;
}

bool WorkSpace::WriteSeeds(std::string filename)
{
	std::ofstream ofile;
	ofile.open(filename);
	current = list_head;
	while (NULL != current)
	{
		PointPair seed = current->GetSeed();
		ofile << seed.x << " " << seed.y << "\n";
		current = current->GetNext();
	}
	ofile.close();
	return true;
}

bool WorkSpace::WriteSuperPixelCSV(std::string filename)
{
	std::ofstream ofile;
	//	std::vector<int> sequence;
	ofile.open(filename);
	current = list_head;
	while (NULL != current)
	{
		ofile << current->SaveState();
		//		ofile << current->write_data();
		current = current->GetNext();
	}
	ofile.close();
	return true;
}

bool WorkSpace::WriteSuperPixels(std::string filename)
{
	std::ofstream ofile;
	std::vector<int> sequence;
	ofile.open(filename, std::ios::binary);
	current = list_head;
	while (NULL != current)
	{
		int identifier = current->GetIdentifier();
		PointPair seed = current->GetSeed();
		RectQuad box = current->GetWindow();
		Color color = current->GetAveColor();

		ofile.write(reinterpret_cast<const char*>(&identifier), sizeof(identifier));
		ofile.write(reinterpret_cast<const char*>(&seed.x), sizeof(seed.x));
		ofile.write(reinterpret_cast<const char*>(&seed.y), sizeof(seed.y));
		ofile.write(reinterpret_cast<const char*>(&box.x0), sizeof(box.x0));
		ofile.write(reinterpret_cast<const char*>(&box.y0), sizeof(box.y0));
		ofile.write(reinterpret_cast<const char*>(&box.x1), sizeof(box.x1));
		ofile.write(reinterpret_cast<const char*>(&box.y1), sizeof(box.y1));
		ofile.write(reinterpret_cast<const char*>(&color.channel[0]), sizeof(color.channel[0]));
		ofile.write(reinterpret_cast<const char*>(&color.channel[1]), sizeof(color.channel[1]));
		ofile.write(reinterpret_cast<const char*>(&color.channel[2]), sizeof(color.channel[2]));
		sequence = current->write_data();
		std::vector<int>::iterator it;
		for (it = sequence.begin(); it != sequence.end(); ++it)
		{
			int num = *it;
			ofile.write(reinterpret_cast<const char*>(&num), sizeof(num));
		}
		sequence.clear();
		current = current->GetNext();
	}
	ofile.close();
	return true;
}

bool WorkSpace::WriteSuperPixelsSVG(std::string filename, int mode, bool polygon, bool fine, int palette_size)
{
	// Mode: 0=normal
	//       1=post-processed
	// Polygon - Use the optimal polygon points.

	bool use_palette = false;
	if (palette_size > 0)
	{
		use_palette = true;
		FindPalette(mode, palette_size);
	}

	std::ofstream ofile;
	//	std::vector<int> sequence;
	ofile.open(filename);
	//ofile << "<svg xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\">";
	ofile << "<?xml version=\"1.0\" standalone=\"no\"?>\n";
	ofile << "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\"\n";
	ofile << "\"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\" >\n";
	ofile << "<svg viewBox = \"0 0 " << (int)(width + 1) << " " << (int)(height + 1) << "\"\n";
	ofile << "xmlns = \"http://www.w3.org/2000/svg\" version = \"1.1\">\n";

	int palette_index = 0;
	if (use_palette)
	{
		palette_index = palette_size - 1;
	}

	while (palette_index >= 0)
	{
		if (0 == mode)
		{
			current = list_head;
		}
		else if (1 == mode)
		{
			current = processed_head;
		}
		else
		{
			return false;  // Should not be used for skeletons. 
		}
		if (use_palette)
		{
			ofile << "<g fill=\"#" << palette->TranslateColorToString(palette->GetColorByIndex(palette_index)) << "\">\n";
		}
		while (NULL != current)
		{
			Path* p = current->GetPathHead();
			while (NULL != p)
			{
				if ((false == use_palette) || (current->GetColorBucket() == palette_index))
				{
					if (use_palette)
					{
						ofile << "<path d = \"M ";
					}
					else {
						ofile << "<path fill=\"#" << p->GetColor() << "\" d = \"M ";
					}
					if (polygon)
					{
						std::vector<Corner> points = p->GetCurve();
						int n = points.size();
						bool first = true;
						Corner p0, p1, p2;  // p1 is the current point, p0 is the previous point, p2 is the following point.
						for (int i = 0; i < n; i++)
						{
							if (i > 0)
							{
								p0 = points[i - 1];
							}
							else {
								p0 = points[n - 1];
							}
							p1 = points[i];
							if (i < (n - 1))
							{
								p2 = points[i + 1];
							}
							else {
								p2 = points[0];
							}
							if (first)
							{
								if (p0.smooth)
								{
									ofile << p1.p0.x << " " << p1.p0.y << " ";  // Start out at first mid point.
								}
								else {
									ofile << p0.c0.x << " " << p0.c0.y << " ";  // Start out at corner before first mid point.
								}
							}

							if (p1.smooth)
							{
								if (false == p0.smooth)
								{
									ofile << "L " << p1.p0.x << " " << p1.p0.y << " ";
								}
								ofile << "C " << p1.c0.x << " " << p1.c0.y << " " << p1.c1.x << " " << p1.c1.y << " " << p1.p1.x << " " << p1.p1.y << " ";
							}
							else {
								ofile << "L " << p1.c0.x << " " << p1.c0.y << " ";
							}
							first = false;
						}
						ofile << "Z\n";
						Path* sub_p = p->GetSubpathHead();
						while (NULL != sub_p)
						{
							ofile << "M ";
							std::vector<Corner> s_points = sub_p->GetCurve();;
							int s_n = s_points.size();
							first = true;
							for (int s_i = 0; s_i < s_n; s_i++)
							{
								if (s_i > 0)
								{
									p0 = s_points[s_i - 1];
								}
								else {
									p0 = s_points[s_n - 1];
								}
								p1 = s_points[s_i];
								if (s_i < (s_n - 1))
								{
									p2 = s_points[s_i + 1];
								}
								else {
									p2 = s_points[0];
								}
								if (first)
								{
									if (p0.smooth)
									{
										ofile << p1.p0.x << " " << p1.p0.y << " ";  // Start out at first mid point.
									}
									else {
										ofile << p0.c0.x << " " << p0.c0.y << " ";  // Start out at corner before first mid point.
									}
								}

								if (p1.smooth)
								{
									if (false == p0.smooth)
									{
										ofile << "L " << p1.p0.x << " " << p1.p0.y << " ";
									}
									ofile << "C " << p1.c0.x << " " << p1.c0.y << " " << p1.c1.x << " " << p1.c1.y << " " << p1.p1.x << " " << p1.p1.y << " ";
								}
								else {
									ofile << "L " << p1.c0.x << " " << p1.c0.y << " ";
								}
								first = false;
							}
							ofile << "Z\n";
							sub_p = sub_p->GetNext();
						}
					}
					else {
						std::vector<PointPair> points;
						if (fine)
						{
							points = p->GetFinePoints();
						}
						else {
							points = p->GetPoints();
						}
						std::vector<PointPair>::iterator p_it;
						bool first = true;
						for (p_it = points.begin(); p_it != points.end(); ++p_it)
						{
							if (false == first)
							{
								ofile << "L ";
							}
							PointPair temp = *p_it;
							ofile << temp.x << " " << temp.y << " ";
							first = false;
						}
						ofile << "Z\n";
						Path* sub_p = p->GetSubpathHead();
						while (NULL != sub_p)
						{
							ofile << "M ";
							std::vector<PointPair> s_points;
							if (fine)
							{
								s_points = sub_p->GetFinePoints();
							}
							else {
								s_points = sub_p->GetPoints();
							}
							std::vector<PointPair>::iterator s_p_it;
							first = true;
							for (s_p_it = s_points.begin(); s_p_it != s_points.end(); ++s_p_it)
							{
								if (false == first)
								{
									ofile << "L ";
								}
								PointPair temp = *s_p_it;
								ofile << temp.x << " " << temp.y << " ";
								first = false;
							}
							ofile << "Z\n";
							sub_p = sub_p->GetNext();
						}
					}
				}
				ofile << "\"\/>\n";
				p = p->GetNext();
			}
			current = current->GetNext();
		}
		if (use_palette)
		{
			ofile << "</g>\n";
		}
		palette_index -= 1;
	}
	ofile << "</svg>";
	ofile.close();
	return true;
}

bool WorkSpace::WritePaintCurvesSVG(std::string filename)
{
	std::ofstream ofile;
	//	std::vector<int> sequence;
	ofile.open(filename);
	//ofile << "<svg xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\">";
	ofile << "<?xml version=\"1.0\" standalone=\"no\"?>\n";
	ofile << "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\"\n";
	ofile << "\"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\" >\n";
	ofile << "<svg viewBox = \"0 0 " << (int)(width + 1) << " " << (int)(height + 1) << "\"\n";
	ofile << "xmlns = \"http://www.w3.org/2000/svg\" version = \"1.1\">\n";

	current = skeleton_head;

	while (NULL != current)
	{
		Path* p = current->GetPathHead();
		while (NULL != p)
		{
			ofile << "<path d = \"M ";
			std::vector<Corner> curve = p->GetCurve();
			int n = curve.size();
			bool first = true;
			Corner p1;  // p1 is the current point, p0 is the previous point, p2 is the following point.
			for (int i = 0; i < n; i++)
			{
				p1 = curve[i];
				if (first)
				{
					ofile << p1.p0.x << " " << p1.p0.y << " ";  // Start out at first mid point.
				}

				if (p1.smooth)
				{
					ofile << "C " << p1.c0.x << " " << p1.c0.y << " " << p1.c1.x << " " << p1.c1.y << " " << p1.p1.x << " " << p1.p1.y << " ";
				}
				else {
					ofile << "L " << p1.c0.x << " " << p1.c0.y << " ";
					ofile << "L " << p1.p1.x << " " << p1.p1.y << " ";
				}
				first = false;
			}
			ofile << "\" ";
			ofile << "fill=\"none\" stroke=\"blue\" stroke-width=\"1\"";
			ofile << "\/>\n";
			p = p->GetNext();
		}

		current = current->GetNext();
	}

	ofile << "</svg>";
	ofile.close();
	return true;
}

bool WorkSpace::SplitSuperPixels(float num_sigmas)
{
	// Split up SuperPixels with large color errors.
	double overall_ave_error = 0;
	current = list_head;
	while (NULL != current)
	{
		overall_ave_error += current->GetAveError();
		current = current->GetNext();
	}
	overall_ave_error = overall_ave_error / list_head->Count();
	double sigma = 0;
	current = list_head;
	while (NULL != current)
	{
		sigma += (current->GetAveError() - overall_ave_error) * (current->GetAveError() - overall_ave_error);
		current = current->GetNext();
	}
	sigma = sqrt(sigma / list_head->Count());
	float limit = overall_ave_error + (num_sigmas * sigma);
	std::cout << "Error: " << overall_ave_error << ", " << sigma << ", " << limit << "\n";
	current = list_head;
	while (NULL != current)
	{
		if (current->GetAveError() > limit)
		{
			PointPair Add_Seed = current->Split(image, diagonals);
			std::cout << current->GetIdentifier() << ", " << Add_Seed.x << ", " << Add_Seed.y << "\n";
			if ((Add_Seed.x != 0) || (Add_Seed.y != 0))
			{
				int Add_Id = 1 + current->GetTail()->GetIdentifier();
				list_tail = new SuperPixel(Add_Id, edge, pixeldata, Add_Seed, NULL, current->GetTail(), this, SPType_Plain);
				InsertSPIntoIndex(list_tail);
			}
		}
		current = current->GetNext();
	}
	return true;
}

bool WorkSpace::ThinSuperPixels(bool glitch3)
{
	SuperPixel* skeleton_current = NULL;
	if (NULL != list_head)
	{
		if (NULL != skeleton_head)
		{
			skeleton_current = skeleton_head;
			while (NULL != skeleton_current)
			{
				if (NULL != skeleton_current->GetNext())
				{
					skeleton_current = skeleton_current->GetNext();
					delete skeleton_current->GetPrev();
				}
				else {
					delete skeleton_current;
					skeleton_current = NULL;
				}
			}
		}
		if (NULL != skeleton_pixeldata)
		{
			delete skeleton_pixeldata;
		}

		if (NULL == processed_head)
		{
			skeleton_pixeldata = new SPixelData(*pixeldata);
			current = list_head;
		}
		else {
			// Because there has already been postprocessing, SuperPixels may be disconnected.  Split them into new SuperPixels, if necessary.
			current = processed_head;
			while (NULL != current)
			{
				if (current->GetPathHead() != current->GetPathTail())
				{
					Path* local_path = current->GetPathHead()->GetNext();
					while (NULL != local_path)
					{
						int point = *(local_path->GetPointSet().begin());
						SeparateSuperPixel(current->GetIdentifier(), point);
						local_path = current->GetPathHead()->GetNext();
					}
				}
				current = current->GetNext();
			}
			skeleton_pixeldata = new SPixelData(*processed_pixeldata);
			current = processed_head;
		}

		skeleton_current = current->Thin(edge, skeleton_pixeldata, NULL, glitch3);
		skeleton_head = skeleton_current;
		while (NULL != current->GetNext())
		{
			current = current->GetNext();
			skeleton_current = current->Thin(edge, skeleton_pixeldata, skeleton_current, glitch3);
		}
		skeleton_tail = skeleton_head->GetTail();
		return true;
	}
	return false;
}

bool WorkSpace::FindPaths(int mode, bool polygon, bool fine)
{
	FindNeighbors(mode);
	if (0 == mode)
	{
		current = list_head;
	}
	else if (1 == mode)
	{
		current = processed_head;
	}
	else
	{
		return false; // Should not be used for skeletons or other things.
	}
	while (NULL != current)
	{
		current->FindPaths(polygon, polygon, fine); // Use meeting points only if creating polygons.
		current = current->GetNext();
	}
	return true;
}

bool WorkSpace::FindNeighbors(int mode)
{
	if (0 == mode)
	{
		current = list_head;
	}
	else if (1 == mode)
	{
		current = processed_head;
	}
	else
	{
		current = skeleton_head;
	}
	while (NULL != current)
	{
		current->SetNeighbors(); // Also finds EdgePixels.
		current = current->GetNext();
	}
	return true;
}

bool WorkSpace::InsertSPIntoIndex(SuperPixel* sp)
{
	SuperPixelIndex index;
	index.identifier = sp->GetIdentifier();
	index.sp = sp;
	superpixel_set.insert(index);
	return true;
}

SuperPixel* WorkSpace::GetByIdentifier(int key)
{
	SuperPixel* ret = NULL;
	SuperPixelIndex index;
	index.identifier = key;
	std::set<SuperPixelIndex>::iterator it;
	it = superpixel_set.find(index);
	if (it != superpixel_set.end())
	{
		index = *it;
		ret = index.sp;
	}
	return ret;
}

bool WorkSpace::RecalcSuperPixelSet()
{
	superpixel_set.clear();
	if (NULL == processed_head)
	{
		current = list_head;
	}
	else {
		current = processed_head;
	}
	while (current != NULL)
	{
		InsertSPIntoIndex(current);
		current = current->GetNext();
	}
	return true;
}

bool WorkSpace::CalculateSizes()
{
	if (NULL == processed_head)
	{
		current = list_head;
	}
	else {
		current = processed_head;
	}

	while (NULL != current)
	{
		current->CalculateSize();
		current = current->GetNext();
	}
	return true;
}

bool WorkSpace::FindPalette(int mode, int palette_size)
{
	// Mode: 0=normal
	//       1=post-processed
	//       2=skeleton

	if (NULL != palette)
	{
		delete(palette);
	}
	if (0 == mode)
	{
		current = list_head;
	}
	else if (1 == mode)
	{
		current = processed_head;
	}
	else
	{
		current = skeleton_head;
	}

	palette = new ColorPalette(palette_size, current);
	if (NULL == palette)
	{
		throw std::runtime_error("Failed to generate ColorPalette.\n");
		return false;
	}
	return true;
}

bool WorkSpace::ReduceToPalette(int mode, int palette_size)
{
	bool ret = FindPalette(mode, palette_size);
	if (false == ret)
	{
		throw std::runtime_error("Unable to reduce to palette.\n");
	}
	if (0 == mode)
	{
		current = list_head;
	}
	else if (1 == mode)
	{
		current = processed_head;
	}
	else
	{
		current = skeleton_head;
	}
	while (NULL != current)
	{
		Color c = current->GetAveColor();
		c = palette->GetColorByIndex(current->GetColorBucket());
		current->SetAveColor(c);
		current = current->GetNext();
	}
	return ret;
}

bool WorkSpace::SeparateSuperPixel(int current_id, int point)
{
	// Step 1: Generate new identifier.

	int new_id = 0;
	if ((NULL == processed_head) || (NULL == processed_pixeldata))
	{
		throw std::runtime_error("SeparateSuperPixel called without processed SuperPixels.\n");
		return false;
	}
	current = processed_head;
	while (NULL != current)
	{
		if (current->GetIdentifier() > new_id)
		{
			new_id = current->GetIdentifier();
		}
		current = current->GetNext();
	}
	new_id++; // Next available number is the new identifier.

	// Step 2: Starting at point, flood-fill the existing processed_pixeldata at point to be new identifier.

	if (false == processed_pixeldata->FloodReplace(point, current_id, new_id))
	{
		throw std::runtime_error("SeparateSuperPixel called on point not in original set.\n");
		return false;
	}

	// Step 3: Calculate new bounding boxes.

	RectQuad orig_box = { 0,0,0,0 };
	RectQuad new_box = { 0,0,0,0 };
	bool orig_seen = false;
	bool new_seen = false;
	for (int j = bbox.y0; j <= bbox.y1; j++)
	{
		for (int i = bbox.x0; i <= bbox.x1; i++)
		{
			int value = processed_pixeldata->GetPixel(i, j);
			if (current_id == value)
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
		throw std::runtime_error("Failed to find new pixels in SeparateSuperPixel.\n");
		return false;
	}
	// Step 4: Confirm that there are original identifier pixels left.

	if (false == orig_seen)
	{
		// Replace original entirely with the new one.
		current = processed_head->GetByIdentifier(current_id);
		if (NULL == current)
		{
			throw std::runtime_error("SeparateSuperPixel called on non-existent SuperPixel.\n");
			return false;
		}
		current->SetIdentifierandBox(new_id, new_box);
	}
	else {
		// Step 5: Create new SuperPixel.

		PointPair seed;
		seed.y = point / width;
		seed.x = point - (seed.y * width);
		SuperPixel* newSP = new SuperPixel(new_id, edge, processed_pixeldata, seed, NULL, processed_tail, this, SPType_Processed);
		if (NULL == newSP)
		{
			throw std::runtime_error("Failed to create new SuperPixel in SeparateSuperPixel.\n");
			return false;
		}
		newSP->SetWindow(new_box);
		InsertSPIntoIndex(newSP);
		processed_tail = newSP;


		// Step 6: Determine which paths (with embedded sub-paths) go with the new SuperPixel.

		// Need to ensure that the processed SuperPixels will be returned, so recompute the SuperPixelSet.
		RecalcSuperPixelSet();
		current = processed_head->GetByIdentifier(current_id);
		Path* local_path = current->GetPathHead();
		while (NULL != local_path)
		{
			std::set<int> points = local_path->GetPointSet();
			int pos = *points.begin();
			int y = pos / width;
			int x = pos - (y * width);
			int value = processed_pixeldata->GetPixel(x, y);
			if (new_id == value) // This path needs to move over to the new SuperPixel.
			{
				local_path->MoveSuperPixel(newSP);
			}
			local_path = local_path->GetNext();
		}

		// Step 7: Set the color of the new SuperPixel to the same as the original.

		newSP->SetAveColor(current->GetAveColor());
	}
	return true;
}

ImageData* WorkSpace::GenerateImage(int mode, Paint_Properties prop)
// Mode: 0=normal
//       1=post-processed
//       2=skeleton
//		 3=painted paths w/out outlines
//       4=painted paths w/ outlines
// Background: 0=white
//             1=black
// glitch1: Create a glitched effect for the straight painted portions.
// glitch2: Use only unsigned char when painting, at low paint flow quantization errors cause values to be lower than they should be.
// sub_pixel: Use more precise calculations to position bristles with sub-pixel accuracy.
{
	Color bg;
	if (0 == prop.background)
	{
		bg.channel[0] = 255;
		bg.channel[1] = 255;
		bg.channel[2] = 255;
	}
	else {
		bg.channel[0] = 0;
		bg.channel[1] = 0;
		bg.channel[2] = 0;
	}
	Color black;
	black.channel[0] = 0;
	black.channel[1] = 0;
	black.channel[2] = 0;
	Color white;
	white.channel[0] = 255;
	white.channel[1] = 255;
	white.channel[2] = 255;

	if (0 == mode)
	{
		return pixeldata->GenerateImage(list_head, bg);
	}
	else if (1 == mode)
	{
		return processed_pixeldata->GenerateImage(processed_head, bg);
	}
	else if (2 == mode)
	{

		current = skeleton_head;
		while (current != NULL)
		{
			std::set<int> points = current->GetSkeletonEndpoints();
			std::set<int>::iterator it;
			for (it = points.begin(); it != points.end(); ++it)
			{
				int pos = *it;
				PointPair xy = current->Pos2XY(pos);
				skeleton_pixeldata->SetPixel(xy.x, xy.y, -1);
			}
			points = current->GetSkeletonIntersections();
			for (it = points.begin(); it != points.end(); ++it)
			{
				int pos = *it;
				PointPair xy = current->Pos2XY(pos);
				skeleton_pixeldata->SetPixel(xy.x, xy.y, -2);
			}
			current = current->GetNext();
		}
		return skeleton_pixeldata->GenerateImage(skeleton_head, bg);
	}
	else if (3 == mode)
	{
		int scale_width = width * prop.paint_scale;
		int scale_height = height * prop.paint_scale;
		ImageData* img = NULL;
		unsigned char* img_data = (unsigned char*)malloc(sizeof(unsigned char) * scale_width * scale_height * 3);
		if (NULL == img_data)
		{
			throw std::runtime_error("Unable to allocate memory for paint image.\n");
			return img;
		}
		img = new ImageData(img_data, scale_width, scale_height, 3, true);
		img->SetBackground(bg);
		// Skeleton path painting
		current = skeleton_head;
		Color second = white;
		SPixelData* local_pixdata = NULL;
		float ave_radius;
		if (NULL != processed_pixeldata)
		{
			local_pixdata = processed_pixeldata;
		}
		else {
			local_pixdata = pixeldata;
		}
		while (current != NULL)
		{
			Path* curve_path = current->GetPathHead();
			while (NULL != curve_path)
			{
				if (false == prop.mix_paints)
				{
					Color lab = img->CIELABconvert(current->GetAveColor());
					lab.channel[0] -= 10.0;
					if (lab.channel[0] > 100.0)
					{
						lab.channel[0] = 100.0;
					}
					else if (lab.channel[0] < 0.0)
					{
						lab.channel[0] = 0.0;
					}
					second = img->RGBconvert(lab);
				}
				std::vector<Corner> curve = curve_path->GetCurve();
				ave_radius = local_pixdata->CalculateRadius(curve.begin(), curve.end(), current->GetIdentifier());
				img->CreateBrush({ 100.0, 100.0 }, current->GetAveColor(), second, ave_radius, prop);
				img->PaintCurve(curve, local_pixdata, current->GetIdentifier(), true);
				if (prop.mix_paints)
				{
					second = current->GetAveColor();
				}
				curve_path = curve_path->GetNext();
			}
			current = current->GetNext();
		}

		if (prop.outline)
		{
			// Outline path painting
			if (NULL != processed_head)
			{
				current = processed_head;
			}
			else {
				current = list_head;
			}
			while (current != NULL)
			{
				Path* curve_path = current->GetPathHead();
				while (NULL != curve_path)
				{
					std::vector<Corner> curve = curve_path->GetCurve();
					Paint_Properties outline_prop = prop;
					outline_prop.flow = 10.0;
					outline_prop.bristles = 10.0;
					outline_prop.flow_variation = 0;
					img->CreateBrush({ 100.0, 100.0 }, black, black, 3, prop);
					img->PaintCurve(curve, NULL, 0, false);
					curve_path = curve_path->GetNext();
				}
				current = current->GetNext();
			}
		}
		img->CollapseWideData();
		return img;
	}
}

ImageData* WorkSpace::Gradient2Image(int mode)
//  Mode:  0=Gradient (edge)
//         1=erode
//         2=dilate
//         3=processed gray
//         4=gray
{
	if ((mode < 0) || (mode > 4))
	{
		return NULL;
	}


	if (0 == mode)
	{
		return edge->Gradient2Image(mode);
	}
	else if (1 == mode)
	{
		return erode->Gradient2Image(mode);
	}
	else if (2 == mode)
	{
		return dilate->Gradient2Image(mode);
	}
	else if (3 == mode)
	{
		return processed_gray->Gradient2Image(mode);
	}

	return gray->Gradient2Image(mode);
}

ColorPalette::ColorPalette(int num, SuperPixel* sp)
{
	MCColorBucket** MCbuckets;
	SuperPixel* current;
	int local_color[3];
	int MC_count;
	unsigned long long color2[3];
	unsigned long long color2_old[3];
	unsigned long long centroid_count = 0;

	size = num;

	MCbuckets = (MCColorBucket**)malloc(sizeof(MCColorBucket*) * size);
	if (NULL == MCbuckets)
	{
		throw std::runtime_error("Unable to allocate memory for MC color buckets.\n");
		return;
	}
	MCbuckets[0] = (MCColorBucket*)malloc(sizeof(MCColorBucket));
	if (NULL == MCbuckets[0])
	{
		throw std::runtime_error("Unable to allocate memory for first MC color bucket.\n");
		free(MCbuckets);
		return;
	}
	// First, set up a single MC bucket representing the full color space.
	MC_count = 0;
	for (int i = 0; i < 3; i++)
	{
		MCbuckets[MC_count]->max[i] = 0;  // Initialize max to zero and min to 255.
		MCbuckets[MC_count]->min[i] = 255;
		color2[i] = 0;
	}

	current = sp->GetHead();
	centroid_count = 0;

	while (NULL != current)
	{
		current->CalculateSize();
		if (current->GetSize() > 0)
		{
			centroid_count++;
			Color c = current->GetAveColor();
			local_color[0] = c.channel[0];
			local_color[1] = c.channel[1];
			local_color[2] = c.channel[2];

			for (int i = 0; i < 3; i++)
			{
				if (local_color[i] > MCbuckets[MC_count]->max[i])
				{
					MCbuckets[MC_count]->max[i] = local_color[i];
				}
				if (local_color[i] < MCbuckets[MC_count]->min[i])
				{
					MCbuckets[MC_count]->min[i] = local_color[i];
				}

				color2[i] += local_color[i] * local_color[i];
			}
			current->SetColorBucket(MC_count);  // Assign SuperPixel to the first MCbucket.
		}
		current = current->GetNext();
	}
	if (centroid_count > 0)
	{
		for (int i = 0; i < 3; i++)
		{
			MCbuckets[MC_count]->centroid[i] = sqrt((double)color2[i] / centroid_count); // Calculate the centroid.
		}
	}
	else {
		size = 0; // No superpixels with any area, so no ColorPalette entries.
		return;
	}
	// Now, create new MCbuckets until we reach the right sizeber.
	while (MC_count < (size - 1))
	{
		MC_count++;
		// Examine existing buckets to find the one with the widest range in any color channel.
		int max_range = 0;  // The largest range.
		int max_MC = 0;     // The MCbucket with the largest range.
		int max_channel = 0; // The channel (0-3, rgb) in which the largest range shows up.
		for (int MCB = 0; MCB < MC_count; MCB++)
		{
			for (int i = 0; i < 3; i++)
			{
				int range = MCbuckets[MCB]->max[i] - MCbuckets[MCB]->min[i];
				if (range > max_range)
				{
					max_range = range;
					max_MC = MCB;
					max_channel = i;
				}
			}
		}

		// Superpixels in the bucket above the centroid value in the selected range stay, other move to the new bucket.
		MCbuckets[MC_count] = (MCColorBucket*)malloc(sizeof(MCColorBucket));
		if (NULL == MCbuckets[MC_count])
		{
			throw std::runtime_error("Unable to allocate memory for MC color bucket.\n");
			for (int i = 0; i < MC_count; i++)
			{
				free(MCbuckets[i]);
			}
			free(MCbuckets);
			return;
		}
		int cutoff = MCbuckets[max_MC]->centroid[max_channel];
		int new_count = 0;
		int old_count = 0;
		for (int i = 0; i < 3; i++)
		{
			MCbuckets[MC_count]->max[i] = 0;  // Initialize max to zero and min to 255.
			MCbuckets[MC_count]->min[i] = 255;
			MCbuckets[max_MC]->max[i] = 0;   // Initialize for the bucket being split.
			MCbuckets[max_MC]->min[i] = 255;
			color2[i] = 0;      // Bucket being created.
			color2_old[i] = 0;  // Bucket being split.
		}
		current = sp->GetHead();
		while (NULL != current)
		{
			if ((current->GetSize() > 0) && (max_MC == current->GetColorBucket()))  // Don't look at non-zero size SuperPixels identified to the right bucket.
			{
				Color c = current->GetAveColor();
				local_color[0] = c.channel[0];
				local_color[1] = c.channel[1];
				local_color[2] = c.channel[2];
				if (local_color[max_channel] < cutoff)  // Test for SuperPixels that need to go to the new bucket.
				{
					current->SetColorBucket(MC_count);
					for (int i = 0; i < 3; i++)
					{
						if (local_color[i] > MCbuckets[MC_count]->max[i])
						{
							MCbuckets[MC_count]->max[i] = local_color[i];
						}
						if (local_color[i] < MCbuckets[MC_count]->min[i])
						{
							MCbuckets[MC_count]->min[i] = local_color[i];
						}
						color2[i] += local_color[i] * local_color[i];
					}
					new_count++;
				}
				else {  // If not going to the new bucket, calculate new centroid numbers for this bucket.
					for (int i = 0; i < 3; i++)
					{
						if (local_color[i] > MCbuckets[max_MC]->max[i])
						{
							MCbuckets[max_MC]->max[i] = local_color[i];
						}
						if (local_color[i] < MCbuckets[max_MC]->min[i])
						{
							MCbuckets[max_MC]->min[i] = local_color[i];
						}
						color2_old[i] += local_color[i] * local_color[i];
					}
					old_count++;
				}
			}
			current = current->GetNext();
		}
		for (int i = 0; i < 3; i++)
		{
			if (new_count > 0)
			{
				MCbuckets[MC_count]->centroid[i] = sqrt((double)color2[i] / new_count); // Calculate centroid for the new bucket.
			}
			if (old_count > 0)
			{
				MCbuckets[max_MC]->centroid[i] = sqrt((double)color2_old[i] / old_count);  // Calculate centroid for old bucket.
			}
		}
	}
	centroid_color.clear();
	for (MC_count = 0; MC_count < size; MC_count++)
	{
		Color c;
		c.channel[0] = MCbuckets[MC_count]->centroid[0];
		c.channel[1] = MCbuckets[MC_count]->centroid[1];
		c.channel[2] = MCbuckets[MC_count]->centroid[2];
		centroid_color.push_back(c);
	}
}

ColorPalette::~ColorPalette()
{
	centroid_color.clear();
}

Color ColorPalette::GetColorByIndex(int index)
{
	if ((index >= 0) && (index < size))
	{
		return centroid_color[index];
	}
	Color c;
	c.channel[0] = 0;
	c.channel[1] = 0;
	c.channel[2] = 0;
	return c;
}

Color ColorPalette::TranslateColor(Color c)  // Takes any color in the colorspace and returns the palette color.
{
	int MC_count = 0;
	int min_dist = 200000;
	int min_MC = 0;
	while (MC_count < size)
	{
		long dist = ((int)c.channel[0] - (int)centroid_color[MC_count].channel[0]) * ((int)c.channel[0] - (int)centroid_color[MC_count].channel[0]);
		dist += ((int)c.channel[1] - (int)centroid_color[MC_count].channel[1]) * ((int)c.channel[1] - (int)centroid_color[MC_count].channel[1]);
		dist += ((int)c.channel[2] - (int)centroid_color[MC_count].channel[2]) * ((int)c.channel[2] - (int)centroid_color[MC_count].channel[2]);
		if (dist < min_dist)
		{
			min_dist = dist;
			min_MC = MC_count;
		}
		MC_count++;
	}
	return centroid_color[min_MC];
}

std::string ColorPalette::TranslateColorToString(Color c)
{
	Color local_color = TranslateColor(c);
	std::stringstream sstream;
	if (local_color.channel[0] < 16)
	{
		sstream << "0";
	}
	sstream << std::hex << (int)local_color.channel[0];
	if (local_color.channel[1] < 16)
	{
		sstream << "0";
	}
	sstream << std::hex << (int)local_color.channel[1];
	if (local_color.channel[2] < 16)
	{
		sstream << "0";
	}
	sstream << std::hex << (int)local_color.channel[2];
	return sstream.str();
}