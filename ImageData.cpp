// Copyright (c) 2023-2025 Steve Young
// Licensed under the MIT License

#include "ImageData.h"
#include "GradData.h"
#include "SPixelData.h"
#include "SuperPixel.h"

#ifdef USE_CUDA
ImageData::ImageData(unsigned char* data_in, int w, int h, int n, bool frac_values, bool watercolor)
#else
ImageData::ImageData(unsigned char* data_in, int w, int h, int n, bool frac_values)
#endif
{
	width = w;
	height = h;
	colorchannels = n;
	brush = NULL;
	// Start with white background color.
	background_color.channel[0] = 255;
	background_color.channel[1] = 255;
	background_color.channel[2] = 255;

	data = (unsigned char*)malloc(sizeof(unsigned char) * width * height * colorchannels);
	if (NULL == data)
	{
		throw (std::runtime_error("Unable to allocate memory for ImageData.\n"));
	}
#ifdef USE_CUDA
	c_device_data = UCharArray(w * colorchannels, h, false);
#endif
	if (frac_values)
	{
		data_wide = (float*)malloc(sizeof(float) * width * height * colorchannels);
		if (NULL == data_wide)
		{
			throw (std::runtime_error("Unable to allocate wide memory for ImageData.\n"));
		}
#ifdef USE_CUDA
		c_device_data_wide = FloatArray(w * colorchannels, h, false);
#endif
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
#ifdef USE_CUDA
		if (!CopyFromHost(data, max_pos, c_device_data))
		{
			throw std::runtime_error("Failed to copy c_device_data to device in ImageData constructor.\n");
		}
#endif
		if (frac_values)
		{
			for (int pos = 0; pos < max_pos; pos++)
			{
				data_wide[pos] = (float)data_in[pos];
			}
#ifdef USE_CUDA
			if (!CopyFromHost(data_wide, max_pos, c_device_data_wide))
			{
				throw std::runtime_error("Failed to copy c_device_data_wide to device in ImageData constructor.\n");
			}
#endif
		}
	}
#ifdef USE_CUDA
	if (watercolor)
	{
		paper = new Paper(width, height, background_color, saturation_dry_value);
	}
	else {
		paper = NULL;
	}
#endif
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
#ifdef USE_CUDA
	if (NULL != paper)
	{
		delete(paper);
	}
	if (NULL != c_device_data)
	{
		FreeUCharArray(c_device_data);
	}
	if (NULL != c_device_data_wide)
	{
		FreeFloatArray(c_device_data_wide);
	}
#endif
}

Color ImageData::GetPixel(int x, int y)
{
	Color ret = { 0, 0, 0 };
	long pos = y * width + x;
	long npos = colorchannels * pos;
	for (int i = 0; (i < 3) && (i < colorchannels); ++i)
	{
		ret.channel[i] = data[npos + i];
	}
	return ret;
}

bool ImageData::CollapseWideData(bool dither)
{
	long pos;
	int line_pixels = width * colorchannels;
	for (long j = 0; j < height; ++j)
	{
		for (long i = 0; i < width; ++i)
		{
			for (int chan = 0; chan < colorchannels; ++chan)
			{
				pos = j * line_pixels + colorchannels * i + chan;
				if (data_wide[pos] >= 254.5)
				{
					data[pos] = 255;
				}
				else if (data_wide[pos] < 0)
				{
					data[pos] = 0;
				}
				else
				{
					data[pos] = (unsigned char)(data_wide[pos] + 0.5);
				}
				if (dither)
				{
					float error = (data_wide[pos] - data[pos]) / 16.0;
					if (i < (width - 1))
					{
						data_wide[pos + colorchannels] += 7.0 * error; // Right pixel
						if (j < (height - 1))
						{
							data_wide[pos + line_pixels + colorchannels] += error; // Lower right pixel
						}
					}
					if (j < (height - 1))
					{
						data_wide[pos + line_pixels] += 5.0 * error;  // Lower pixel
						if (i > 0)
						{
							data_wide[pos + line_pixels - colorchannels] += 3.0 * error; // Lower left pixel
						}
					}
				}
			}
		}
	}
	return true;
}

bool ImageData::CreateBrush(FloatPointPair start, Color c, Color sec, int r, Paint_Properties prop, int pigment_index)
{
	bool ret = true;
	if (NULL != brush)
	{
		delete(brush);
	}
#ifdef USE_CUDA
	brush = new Brush(start, c, sec, r, 10, prop, paper, pigment_index);
#else
	brush = new Brush(start, c, sec, r, 10, prop);
#endif
	if (NULL == brush)
	{
		throw std::runtime_error("Failed to create Brush in ImageData::CreateCudaBrush function.\n");
		ret = false;
	}
	return ret;
}

bool ImageData::PaintCurve(std::vector<Corner> curve, SPixelData* mask, int mask_value, bool use_mask, SPixelData* extinguish_mask)
{
	if (NULL == data_wide)
	{
		return false;
	}
	Corner local_corner;
	std::vector<Corner>::iterator it;
#ifdef USE_CUDA
	if (brush->GetWatercolor())
	{
		brush->StrokeBegin();
	}
#endif
	for (it = curve.begin(); it != curve.end(); ++it)
	{
		local_corner = *it;
		if (local_corner.smooth)
		{
			brush->PaintCorner(local_corner, data_wide, width, height, mask, mask_value, use_mask, extinguish_mask);
		}
		else {
			FloatPointPair dir = { 0, 0 };
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
			brush->PaintTo(local_corner.c0, dir, data_wide, width, height, mask, mask_value, use_mask, local_corner.radius_p0, local_corner.radius_c0, true, extinguish_mask);
			brush->PaintTo(local_corner.p1, dir, data_wide, width, height, mask, mask_value, use_mask, local_corner.radius_c0, local_corner.radius_p1, false, extinguish_mask);
		}
	}
#ifdef USE_CUDA
	if (brush->GetWatercolor())
	{
		brush->StrokeEnd();
	}
#endif
	return true;
}

unsigned char* ImageData::GetData()
{
	return data;
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
#ifdef USE_CUDA
	if (!ResetUCharArray(c_device_data, width, height * colorchannels, 0))
	{
		throw std::runtime_error("Unable to reset c_device_data in ImageData::Reset().\n");
		return false;
	}
#endif
	if (NULL != data_wide)
	{
		for (int pos = 0; pos < max_pos; pos++)
		{
			data_wide[pos] = 0;
		}
#ifdef USE_CUDA
		if (!ResetFloatArray(c_device_data_wide, width, height * colorchannels, 0))
		{
			throw std::runtime_error("Unable to reset c_device_data_wide in ImageData::Reset().\n");
			return false;
		}
#endif
	}
	return true;
}

bool ImageData::SetBackground(Color c)
{
	int max_pos = width * height;
	background_color = c;
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
#ifdef USE_CUDA
	if (NULL != paper)
	{
		paper->SetBackgroundColor(background_color);
	}
#endif
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

#ifdef USE_CUDA
bool ImageData::ProcessWatercolor()
{
	bool ret = true;
	if (NULL == paper)
	{
		throw std::runtime_error("ProcessWatercolor called without pointer to Paper in ImageData.\n");
		ret = false;
	}
	else {
		ret = paper->Process(false);
	}
	return ret;
}

bool ImageData::RenderWatercolor()
{
	bool ret = true;
	if (NULL == paper)
	{
		throw std::runtime_error("Called RenderWaterColor with no point to Paper in ImageData.\n");
		ret = false;
	}
	else {
		paper->Render(data);
	}
	return ret;
}

Paper* ImageData::GetPaper()
{
	return paper;
}
#endif

GradData* ImageData::gen_gray(int channel, int nchannel)
{
	unsigned char* gray;
	GradData* ret = NULL;
	int w = GetWidth();
	int h = GetHeight();
	int n = GetColorChannels();
#ifdef USE_CUDA
	unsigned char* c_device_gray = UCharArray(w, h, false);
#endif
	gray = (unsigned char*)malloc(sizeof(unsigned char) * w * h);
	if (NULL == gray)
	{
		throw (std::runtime_error("Unable to allocate memory for gray image.\n"));
		return NULL;
	}
#ifdef USE_CUDA
	if (!c_gen_gray(c_device_data, w, h, n, channel, nchannel, c_device_gray))
	{
		throw std::runtime_error("Error processing c_gen_gray in ImageData::gen_gray.\n");
		return NULL;
	}
	if (!CopyToHost(c_device_gray, w * h, gray))
	{
		throw std::runtime_error("Error copying data to host in ImageData::gen_gray.\n");
		return NULL;
	}
#else
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
				gray[pos] = c.channel[0];
			}
		}
	}
#endif

	ret = new GradData(gray, w, h);
	if (NULL == ret)
	{
		throw std::runtime_error("Failed to create GradData object in gen_gray.\n");
	}
	free(gray);
#ifdef USE_CUDA
	FreeUCharArray(c_device_gray);
#endif
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


Color ImageData::CIELABconvert(float r, float g, float b)
{
	Color ret = { 0, 0, 0 }; // L*, a*, b* format.
	double lin[3] = { 0, 0, 0 };
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
	double XYZ[3] = { 0, 0, 0 };
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
	Color ret = { 0, 0, 0 };
	double XYZ[3] = { 0, 0, 0 };
	// Calculate Y
	double fXYZ[3] = { 0, 0, 0 };
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
		if (ret.channel[i] <= 0.0)
		{
			ret.channel[i] = 0.0;
		}
		else if (ret.channel[i] >= 1.0)
		{
			ret.channel[i] = 255.0;
		}
		else {
			ret.channel[i] = 255.0 * ret.channel[i];
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

void write_png_to_mem(void* context, void* data, int size)
{
	SuperPixel* sp = (SuperPixel*)context;

	std::string fill_string;
	fill_string.clear();
	int count = 0;
	int remaining_bytes = 0;
	unsigned char input[3] = { 0, 0, 0 };
	unsigned char* input_data = (unsigned char*)data;
	unsigned char output[4] = { 0, 0, 0, 0 };
	while (count < size)
	{
		for (int i = 0; i < 3; ++i)
		{
			if ((count + i) < size)
			{
				input[i] = input_data[i + count];
			}
			else {
				input[i] = 0;
				if (0 == remaining_bytes)
				{
					remaining_bytes = 3 - i;
				}
			}
		}
		uint32_t input_concat = 0;
		for (int i = 0; i < 3; ++i)
		{
			input_concat = input_concat | (input[i] << (16 - 8 * i));
		}
		for (int i = 0; i < 4; ++i)
		{
			output[i] = 63 & (input_concat >> (18 - 6 * i));
		}

		for (int i = 0; i < 4; ++i)
		{
			std::string c;
			if ((remaining_bytes > 0) && (i > 1) && (0 == output[i]) && (i + remaining_bytes >= 4))  // =
			{
				c = (char)61;
			}
			else if (output[i] < 26) // Capital letters
			{
				c = (char)(output[i] + 65);
			}
			else if (output[i] < 52) // Lower case letters
			{
				c = (char)(output[i] + 71);
			}
			else if (output[i] < 62) // Numerals
			{
				c = (char)(output[i] - 4);
			}
			else if (62 == output[i]) // +
			{
				c = (char)43;
			}
			else                      // /
			{
				c = (char)47;
			}
			fill_string.append(c);
		}
		count = count + 3;
	}
	sp->SetFillImage(fill_string);
}
