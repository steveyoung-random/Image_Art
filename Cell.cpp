// Copyright (c) 2023-2025 Steve Young
// Licensed under the MIT License

#include "Cell.h"

Cell::Cell(GradData* graddat, SPixelData* pixdat, int cx, int cy, int w, int h, int buf)
{
	gradientdata = graddat;
	pixeldata = pixdat;
	if ((buf >= h / 2) || (buf >= w / 2))
	{
		buf = 0;
	}
	if ((cx < 0) || (cy < 0) || (cx >= gradientdata->GetWidth()) || (cy >= gradientdata->GetHeight()))
	{
		throw std::runtime_error("Values out of bounds for Cell.\n");
	}
	buffer = buf;
	center_x = cx;
	center_y = cy;
	box.x0 = cx - w / 2;
	box.y0 = cy - h / 2;
	box.x1 = box.x0 + w - 1;
	box.y1 = box.y0 + h - 1;
	if (box.x0 < 0)
	{
		box.x0 = 0;
	}
	if (box.y0 < 0)
	{
		box.y0 = 0;
	}
	if (box.x1 >= gradientdata->GetWidth())
	{
		box.x1 = gradientdata->GetWidth() - 1;
	}
	if (box.y1 >= gradientdata->GetHeight())
	{
		box.y1 = gradientdata->GetHeight() - 1;
	}
	window.x0 = box.x0 + buffer;
	window.y0 = box.y0 + buffer;
	window.x1 = box.x1 - buffer;
	window.y1 = box.y1 - buffer;
	numPixels = 0;
	head = NULL;
	tail = NULL;
}

Cell::~Cell()
{

}

bool Cell::FindSeed(SPixelData* mask, int mask_value, bool diagonals)
{
	bool ignore_mask;
	unsigned char min = 255;
	int count = 0;
	SuperPixel* current = NULL;
	head = NULL;
	tail = NULL;

	ignore_mask = (mask == NULL);

	for (int j = window.y0; j <= window.y1; j++)
	{
		for (int i = window.x0; i <= window.x1; i++)
		{
			if (ignore_mask || (mask_value == mask->GetPixel(i, j)))
			{
				unsigned char pix = gradientdata->GetPixel(i, j);
				if (pix < min)
					min = pix;
			}
		}
	}
	pixeldata->Reset();
	for (int j = window.y0; j <= window.y1; j++)
	{
		for (int i = window.x0; i <= window.x1; i++)
		{
			if ((ignore_mask || (mask_value == mask->GetPixel(i, j)))&&(0 == pixeldata->GetPixel(i, j)))
			{
				if (min == gradientdata->GetPixel(i, j))
				{
					count++;
					PointPair ij;
					ij.x = i;
					ij.y = j;
					if (NULL == head)
					{
						head = new SuperPixel(count, gradientdata, pixeldata, ij, NULL, NULL);
						tail = head;
					}
					else {
						tail = new SuperPixel(count, gradientdata, pixeldata, ij, NULL, tail);
					}
					while (count == tail->Grow(min, true, true, box, mask, mask_value, diagonals));
				}
			}
		}
	}
	if (NULL == head)
	{
		return false;
	}
	while (head != tail)
	{
		min++;
		current = head;
		while (NULL != current)
		{
			current->SetPrevSize();
			current = current->GetNext();
		}
		current = head;
		while (NULL != current)
		{
			bool cont_flag = false;
			int grow_ret = current->Grow(min, true, false, box, mask, mask_value, diagonals);
			while (grow_ret == current->GetIdentifier())
			{
				grow_ret = current->Grow(min, true, false, box, mask, mask_value, diagonals);
			}
			if (grow_ret > 0) // Two SuperPixels have touched at min.  Need to compare sizes at previous step.
			{
				SuperPixel* intersect = head;
				while (intersect->GetIdentifier() != grow_ret)
				{
					intersect = intersect->GetNext();
					if (NULL == intersect)
					{
						throw std::runtime_error("Intersecting SuperPixel not found.\n");
					}
				}
				if (current->GetPrevSize() >= intersect->GetPrevSize())
				{
					current->Absorb(intersect);
				}
				else {
					intersect->Absorb(current);
					current = intersect->GetHead();
				}
				cont_flag = true;
				head = current->GetHead();
				tail = current->GetTail();
			}
			if (false == cont_flag)
			{
				current = current->GetNext();
			}
		}
	}
	if (NULL == head)
	{
		throw std::runtime_error("Failed to find seed for Cell.\n");
	}
	seed = head->GetSeed();
	delete head;
	return true;
}

PointPair Cell::GetSeed()
{
	return seed;
}
