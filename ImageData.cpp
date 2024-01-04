// Copyright (c) 2023-2024 Steve Young
// Licensed under the MIT License

#include "ImageData.h"
#include "GradData.h"
#include "SPixelData.h"

ImageData::ImageData(unsigned char* data_in, int w, int h, int n, bool frac_values)
{
	width = w;
	height = h;
	colorchannels = n;
	data = (unsigned char*)malloc(sizeof(unsigned char) * width * height * colorchannels);
	if (NULL == data)
	{
		throw (std::runtime_error("Unable to allocate memory for ImageData.\n"));
	}
	if (frac_values)
	{
		data_wide = (float*)malloc(sizeof(float) * width * height * colorchannels);
		if (NULL == data_wide)
		{
			throw (std::runtime_error("Unable to allocate wide memory for ImageData.\n"));
		}
	}
	else {
		data_wide = NULL;
	}
	if (data_in == NULL)
	{
		Reset();
	}
	else {
		int max_pos = width * height * colorchannels;
		for (int pos = 0; pos < max_pos; pos++)
		{
			data[pos] = data_in[pos];
		}
		if (frac_values)
		{
			for (int pos = 0; pos < max_pos; pos++)
			{
				data_wide[pos] = (float)data_in[pos];
			}
		}
	}
}

ImageData::~ImageData()
{
	if (NULL != data)
	{
		free(data);
	}
	if (NULL != data_wide)
	{
		free(data_wide);
	}
	if (NULL != brush)
	{
		delete(brush);
	}
}

Color ImageData::GetPixel(int x, int y)
{
	Color ret;
	long pos = y * width + x;
	long npos = colorchannels * pos;
	for (int i = 0; (i < 3) && (i < colorchannels); ++i)
	{
		ret.channel[i] = data[npos + i];
	}
	return ret;
}

bool ImageData::CollapseWideData()
{
	long maxpos = width * height * colorchannels;
	for (long pos = 0; pos < maxpos; ++pos)
	{
		if (data_wide[pos] > 255)
		{
			data[pos] = 255;
		}
		else if (data_wide[pos] < 0)
		{
			data[pos] = 0;
		}
		else
			data[pos] = (unsigned char)data_wide[pos];
	}
	return true;
}

bool ImageData::CreateBrush(FloatPointPair start, Color c, Color sec, int r, Paint_Properties prop)
{
	if (NULL != brush)
	{
		delete(brush);
	}
	brush = new Brush(start, c, sec, r, 10, prop);
	return true;
}

bool ImageData::PaintCurve(std::vector<Corner> curve, SPixelData* mask, int mask_value)
{
	if (NULL == data_wide)
	{
		return false;
	}
	Corner local_corner;
	std::vector<Corner>::iterator it;
	for (it = curve.begin(); it != curve.end(); ++it)
	{
		local_corner = *it;
		if (local_corner.smooth)
		{
			brush->PaintCorner2(local_corner, data_wide, width, height, mask, mask_value);
		}
		else {
			FloatPointPair dir;
			if (curve.begin() == it)
			{
				if ((abs(local_corner.c0.x - local_corner.p0.x) > EFFECTIVE_ZERO) || (abs(local_corner.c0.y - local_corner.p0.y) > EFFECTIVE_ZERO))
				{
					dir.x = local_corner.c0.x - local_corner.p0.x;
					dir.y = local_corner.c0.y - local_corner.p0.y;
				}
				else {
					dir.x = local_corner.p1.x - local_corner.p0.x;
					dir.y = local_corner.p1.y - local_corner.p0.y;
				}
				brush->SetOrientation(dir);
			}
			float o = brush->GetOrientation();
			dir.x = cos(o);
			dir.y = sin(o);
			brush->MoveTo(local_corner.p0);
			brush->PaintTo2(local_corner.c0, dir, data_wide, width, height, mask, mask_value, local_corner.radius_p0, local_corner.radius_c0, true);
			brush->PaintTo2(local_corner.p1, dir, data_wide, width, height, mask, mask_value, local_corner.radius_c0, local_corner.radius_p1, false);
		}
	}
	return true;
}

int ImageData::GetWidth()
{
	return width;
}

int ImageData::GetHeight()
{
	return height;
}

int ImageData::GetColorChannels()
{
	return colorchannels;
}

bool ImageData::SetPixel(int x, int y, Color c)
{
	long pos = y * width + x;
	long npos = colorchannels * pos;
	if ((x >= 0) && (y >= 0) && (x < width) && (y < height))
	{
		for (int i = 0; (i < 3) && (i < colorchannels); ++i)
		{
			data[npos + i] = c.channel[i];
		}
		if (NULL != data_wide)
		{
			for (int i = 0; (i < 3) && (i < colorchannels); ++i)
			{
				data_wide[npos + i] = c.channel[i];
			}
		}
		return true;
	}
	else {
		throw std::runtime_error("Attempt to set pixel outside of bounds.\n");
		return false;
	}
}

bool ImageData::Reset()
{
	int max_pos = width * height * colorchannels;
	for (int pos = 0; pos < max_pos; pos++)
	{
		data[pos] = 0;
	}
	if (NULL != data_wide)
	{
		for (int pos = 0; pos < max_pos; pos++)
		{
			data_wide[pos] = 0;
		}
	}
	return true;
}

bool ImageData::SetBackground(Color c)
{
	int max_pos = width * height;
	for (int pos = 0; pos < max_pos; pos++)
	{
		for (int i = 0; i < colorchannels; i++)
		{
			data[pos * colorchannels + i] = c.channel[i];
		}
		if (NULL != data_wide)
		{
			for (int i = 0; i < colorchannels; i++)
			{
				data_wide[pos * colorchannels + i] = c.channel[i];
			}
		}
	}
	return true;
}

bool ImageData::write_file(std::string filename)
{
	std::stringstream ss;
	ss << filename << ".png";
	std::string s = ss.str();

	if (0 == stbi_write_png(s.c_str(), width, height, colorchannels, data, width * colorchannels))
	{
		throw std::runtime_error("Unable to write out image.\n");
	}
	return true;
}

GradData* ImageData::gen_gray(int channel, int nchannel)
{
	unsigned char* gray;
	GradData* ret = NULL;
	int w = GetWidth();
	int h = GetHeight();
	int n = GetColorChannels();
	gray = (unsigned char*)malloc(sizeof(unsigned char) * w * h);
	if (NULL == gray)
	{
		throw (std::runtime_error("Unable to allocate memory for gray image.\n"));
		return NULL;
	}

	if (n > 2)
	{
		if (0 == channel)
		{
			// Calculate grayscale.
			for (int j = 0; j < h; j++)
			{
				for (int i = 0; i < w; i++)
				{
					Color c = GetPixel(i, j);
					long pos = j * w + i;
					gray[pos] = (unsigned char)(sqrt(c.channel[0] * c.channel[0] + c.channel[1] * c.channel[1] + c.channel[2] * c.channel[2]) * 0.57735);
				}
			}
		}
		else if ((0 < channel) && (channel <= 3) && (0 == nchannel)) {
			for (int j = 0; j < h; j++)
			{
				for (int i = 0; i < w; i++)
				{
					Color c = GetPixel(i, j);
					long pos = j * w + i;
					gray[pos] = c.channel[channel - 1];
				}
			}
		}
		else if ((0 < channel) && (channel <= 3) && (0 < nchannel) && (nchannel <= 3) && (channel != nchannel)) // channel and nchannel both nonzero and in range.
		{
			for (int j = 0; j < h; j++)
			{
				for (int i = 0; i < w; i++)
				{
					Color c = GetPixel(i, j);
					long pos = j * w + i;
					int combined;
					combined = c.channel[channel - 1];
					combined = combined - c.channel[nchannel - 1];
					combined = (combined + 255) / 2;
					gray[pos] = (unsigned char)combined;
				}
			}
		}
		else {
			throw std::runtime_error("Values for channel and nchannel out of range or conflict.\n");
			return ret;
		}
	}
	else {
		// Calculate grayscale.

		for (int j = 0; j < h; j++)
		{
			for (int i = 0; i < w; i++)
			{
				Color c = GetPixel(i, j);
				long pos = j * w + i;
				gray[pos] = c.channel[i];
			}
		}
	}

	ret = new GradData(gray, w, h);
	if (NULL == ret)
	{
		throw std::runtime_error("Failed to create GradData object in gen_gray.\n");
	}
	free(gray);
	return ret;
}

GradData* ImageData::gen_diff(ImageData* image2)
{
	int x = GetWidth();
	int y = GetHeight();
	int n = GetColorChannels();
	unsigned char* diff;
	GradData* ret = NULL;
	diff = (unsigned char*)malloc(sizeof(unsigned char) * x * y);
	if (NULL == diff)
	{
		throw (std::runtime_error("Unable to allocate memory for diff image.\n"));
		return NULL;
	}
	if (n > 2)
	{
		// Calculate diffs.
		for (int j = 0; j < y; j++)
		{
			for (int i = 0; i < x; i++)
			{
				int pos = i + j * x;
				Color c1 = GetPixel(i, j);
				Color c2 = image2->GetPixel(i, j);
				int diff_1 = c1.channel[0] - c2.channel[0];
				int diff_2 = c1.channel[1] - c2.channel[1];
				int diff_3 = c1.channel[2] - c2.channel[2];
				diff[pos] = (unsigned char)(sqrt(diff_1 * diff_1 + diff_2 * diff_2 + diff_3 * diff_3) * 0.57735);
			}
		}
	}
	else {
		// add
		throw std::runtime_error("Unsupported bit depth.\n");
		return NULL;
	}

	ret = new GradData(diff, x, y);
	if (NULL == ret)
	{
		throw std::runtime_error("Failed to create GradData object in gen_diff.\n");
	}
	free(diff);
	return ret;
}


Color ImageData::CIELABconvert(unsigned char r, unsigned char g, unsigned char b)
{
	Color ret; // L*, a*, b* format.
	double lin[3];
	// Step 1, linearize Display P3 values.
	lin[0] = r / 255.0;
	lin[1] = g / 255.0;
	lin[2] = b / 255.0;
	for (int i = 0; i < 3; ++i)
	{
		if (lin[i] <= 0.04045)
		{
			lin[i] = lin[i] / 12.92;
		}
		else {
			lin[i] = pow((lin[i] + 0.055) / 1.055, 2.4);
		}
	}
	// Step 2, convert linearized RGB to XYZ
	double XYZ[3];
	//XYZ[0] = 0.4865709 * lin[0] + 0.2656677 * lin[1] + 0.1982173 * lin[2];
	//XYZ[1] = 0.2289746 * lin[0] + 0.6917385 * lin[1] + 0.0792869 * lin[2];
	//XYZ[2] = 0.0000000 * lin[0] + 0.0451134 * lin[1] + 1.0439444 * lin[2];

	XYZ[0] = 0.4124564 * lin[0] + 0.3575761 * lin[1] + 0.1804375 * lin[2];
	XYZ[1] = 0.2126729 * lin[0] + 0.7151522 * lin[1] + 0.0721750 * lin[2];
	XYZ[2] = 0.0193339 * lin[0] + 0.1191920 * lin[1] + 0.9503041 * lin[2];

	// Step 3, convert XYZ to CIELAB

	XYZ[0] = XYZ[0] / 0.95047;
	//	XYZ[1] = XYZ[1] / 1.000;
	XYZ[2] = XYZ[2] / 1.08883;
	for (int i = 0; i < 3; ++i)
	{
		if (XYZ[i] > 0.008856)
		{
			XYZ[i] = pow(XYZ[i], (1.0 / 3.0));
		}
		else {
			XYZ[i] = 7.787 * XYZ[i] + (16.0 / 116.0);
		}
	}
	ret.channel[0] = 116.0 * XYZ[1] - 16.0;
	ret.channel[1] = 500.0 * (XYZ[0] - XYZ[1]);
	ret.channel[2] = 200.0 * (XYZ[1] - XYZ[2]);
	return ret;
}

Color ImageData::RGBconvert(float L, float a, float b)
{
	Color ret;
	double XYZ[3];
	// Calculate Y
	double fXYZ[3];
	fXYZ[1] = (L + 16.0) / 116.0;
	fXYZ[0] = a / 500.0 + fXYZ[1];
	fXYZ[2] = fXYZ[1] - b / 200.0;
	if (L > 8.856)
	{
		XYZ[1] = pow(fXYZ[1], 3.0);
	}
	else {
		XYZ[1] = L / 903.3;
	}
	// Calculate X
	if (fXYZ[0] > 0.2069)
	{
		XYZ[0] = 0.95047 * pow(fXYZ[0], 3.0);
	}
	else {
		XYZ[0] = 0.95047 * (fXYZ[0] - (4.0 / 29.0)) / 7.787;
	}
	// Calculate Z
	if (fXYZ[2] > 0.2069)
	{
		XYZ[2] = 1.08883 * pow(fXYZ[2], 3.0);
	}
	else {
		XYZ[2] = 1.08883 * (fXYZ[2] - (4.0 / 29.0)) / 7.787;
	}
	// Calculate linearized RGB
	//ret.channel[0] = 2.4934973 * XYZ[0] - 0.9313838 * XYZ[1] - 0.4027109 * XYZ[2];
	//ret.channel[1] = -0.8294893 * XYZ[0] + 1.7626642 * XYZ[1] + 0.0236248 * XYZ[2];
	//ret.channel[2] = 0.0358459 * XYZ[0] - 0.0761724 * XYZ[1] + 0.9568845 * XYZ[2];

	ret.channel[0] = 3.2404548 * XYZ[0] - 1.5371389 * XYZ[1] - 0.4985315 * XYZ[2];
	ret.channel[1] = -0.9692664 * XYZ[0] + 1.8760109 * XYZ[1] + 0.0415561 * XYZ[2];
	ret.channel[2] = 0.0556434 * XYZ[0] - 0.2040259 * XYZ[1] + 1.0572252 * XYZ[2];
	// Gamma correction and 8 bit conversion.
	for (int i = 0; i < 3; ++i)
	{
		if (ret.channel[i] < 0.0031308)
		{
			ret.channel[i] = ret.channel[i] * 12.92;
		}
		else {
			ret.channel[i] = 1.055 * pow(ret.channel[i], (1 / 2.4)) - 0.055;
		}
		if (ret.channel[i] <= 0)
		{
			ret.channel[i] = 0;
		}
		else if (ret.channel[i] >= 1.0)
		{
			ret.channel[i] = 255;
		}
		else {
			ret.channel[i] = 255 * ret.channel[i];
		}
	}
	return ret;
}

Color ImageData::CIELABconvert(Color input)
{
	return CIELABconvert(input.channel[0], input.channel[1], input.channel[2]);
}

Color ImageData::RGBconvert(Color input)
{
	return RGBconvert(input.channel[0], input.channel[1], input.channel[2]);
}