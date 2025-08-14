#include <iostream>
#include <sstream>
#include <string>

#include "Brush.h"
#include "Paper.h"
#include <random>

//#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
//#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#define WIDTH 800
#define HEIGHT 800
#define RADIUS 20
#define BRISTLES 500
#define FLOW 7
#define FLOW_VARIATION 5
#define FILENAME "test.png"
#define PAPERNAME "paper.png"
#define WCOLORTEST "wcolor_test.png"



int main()
{
	int width = WIDTH;
	int height = HEIGHT;
	unsigned char* data;
	Paint_Properties prop;

	try
	{
		std::stringstream ss;
		ss << FILENAME;
		std::string s = ss.str();

		std::stringstream ss_paper;
		ss_paper << PAPERNAME;
		std::string s_paper = ss_paper.str();

		FloatPointPair start = { WIDTH / 8, HEIGHT / 3 };
		FloatPointPair position = { WIDTH * 7 / 8, HEIGHT * 2 / 3 };
		Color color = { 255, 0, 0 };
		Color second = { 0, 255, 0 };
		Color third = { 0, 0, 255 };


		Color paper_color;
		paper_color.channel[0] = 255;
		paper_color.channel[1] = 255;
		paper_color.channel[2] = 255;

		Paper* p = new Paper(WIDTH, HEIGHT, paper_color, saturation_dry_value);
		if (NULL == p)
		{
			throw std::runtime_error("Failed to create Paper instance.\n");
			return -1;
		}

		data = (unsigned char*)malloc(sizeof(unsigned char) * width * height * 3);
		if (NULL == data)
		{
			throw std::runtime_error("Failed to allocate data for image.\n");
			return -1;
		}

		float* thickness = NULL;
		thickness = (float*)malloc(WIDTH * HEIGHT * sizeof(float));
		if (NULL == thickness)
		{
			throw std::runtime_error("Failed to allocate memory for thickness.\n");
			exit(1);
		}
		if (p->GetThickness(thickness))
		{
			for (int j = 0; j < HEIGHT; ++j)
			{
				for (int i = 0; i < WIDTH; ++i)
				{
					long pos = 3 * (i + j * WIDTH);
					data[pos] = (unsigned char)255 * thickness[i + j * WIDTH];
					data[pos + 1] = data[pos];
					data[pos + 2] = data[pos];
				}
			}
			if (0 == stbi_write_png(s_paper.c_str(), width, height, 3, data, width * 3))
			{
				throw std::runtime_error("Unable to write out test image.\n");
			}
		}
		else {
			throw std::runtime_error("Failed to get thickness from Paper.\n");
			exit(1);
		}

		

		//	Pigment* pgmnt = new Pigment(WIDTH, HEIGHT, 1.62, 0.61, 1.64, 0.01, 0.012, 0.003, 0.09, 1.0, 0.41);
		Pigment* pgmnt = new Pigment(WIDTH, HEIGHT, 0.86, 0.86, 0.06, 0.005, 0.005, 0.09, 0.09, 1.0, 0.41);

		if (NULL == pgmnt)
		{
			throw std::runtime_error("Failed to create pigment.\n");
			return -1;
		}
		p->SetPigment(pgmnt);
		for (int count = 0; count < 0.75 * width; ++count)
		{
			p->Dab(width / 8 + count, height / 8 + count, 40, 0.0, 0.02, 0);
			p->SetVelocity(width / 8 + count, height / 8 + count, 1.0, 1.0, true, false);
		}
		p->SetVelocity(0, 0, 0, 0, true, true);
		std::cout << "Brush up. ";
		try {
			p->Process(false);
		}
		catch (std::runtime_error e)
		{
			std::cout << "Exception: " << e.what() << "\n";
			exit(1);
		}
		catch (...)
		{
			std::cout << "Unhandled exception.\n";
			exit(1);
		}
		std::cout << "Done with first pass.\n";
		pgmnt = new Pigment(WIDTH, HEIGHT, 1.52, 0.32, 0.25, 0.06, 0.26, 0.4, 0.01, 1.0, 0.31);
		if (NULL == pgmnt)
		{
			throw std::runtime_error("Failed to create pigment.\n");
			return -1;
		}
		p->SetPigment(pgmnt);
		pgmnt = new Pigment(WIDTH, HEIGHT, 0.10, 0.26, 3.45, 0.97, 0.65, 0.007, 0.05, 6.4, 0.91);
		if (NULL == pgmnt)
		{
			throw std::runtime_error("Failed to create pigment.\n");
			return -1;
		}
		p->SetPigment(pgmnt);
		for (int count = 0; count < 0.75 * width; ++count)
		{
			float frac = (float)count / (float)width;
			p->Dab(width * 7 / 8 - count, height / 8 + count, 40, 0.8, 0.0025 + 0.05 * frac, 1);
			p->SetVelocity(width * 7 / 8 - count, height / 8 + count, -1.0, 1.0, true, false);
		}
		p->SetVelocity(0, 0, 0, 0, true, true);
		std::cout << "Brush up. ";
		std::cout << "Done with second pass.\n";



		for (int count = 0; count < 0.75 * height; ++count)
		{
			p->Dab(width * 3 / 8, height / 8 + count, 40, 0.8, 0.025 + 0.01 * (count / height), 2);
			p->SetVelocity(width * 3 / 8, height / 8 + count, 0.0, 1.0, true, false);
		}
		p->SetVelocity(0, 0, 0, 0, true, true);
		std::cout << "Brush up. ";
		try
		{
			p->Process(false);
		}
		catch (std::runtime_error e)
		{
			std::cout << "Exception: " << e.what() << "\n";
			exit(1);
		}
		catch (...)
		{
			std::cout << "Unhandled exception.\n";
			exit(1);
		}
		// Change colors.
		//p->GetPigments()[0]->UpdateColor(&color, 0.5);
		//p->GetPigments()[1]->UpdateColor(&second, 0.5);
		//p->GetPigments()[2]->UpdateColor(&third, 0.5);

		//int n = 5;
		//int spacing = width / (2 * n + 1);
		//int p_num = 2;
		//std::random_device rd;
		//std::mt19937 gen(rd());

		//std::uniform_real_distribution<> radius(spacing/2, spacing * 2);
		//std::uniform_real_distribution<> move_x(-spacing / 4, spacing / 4);
		//std::uniform_real_distribution<> move_y(-spacing / 4, spacing / 4);
		//std::uniform_int_distribution<> c_val(0, 255);

		//for (int i = 0; i < n; ++i)
		//{
		//	for (int j = 0; j < n; ++j)
		//	{
		//		for (int color_chan = 0; color_chan < 3; ++color_chan)
		//		{
		//			color.channel[color_chan] = c_val(gen);
		//		}
		//		p_num++;
		//		pgmnt = new Pigment(WIDTH, HEIGHT, &color, 0.5, 0.01, 1.0, 0.31);
		//		if (NULL == pgmnt)
		//		{
		//			throw std::runtime_error("Failed to create pigment.\n");
		//			return -1;
		//		}
		//		p->SetPigment(pgmnt);
		//		int x0 = (2 * (i + 1)) * spacing;
		//		int x1 = x0 + move_x(gen);
		//		int y0 = (2 * (j + 1)) * spacing;
		//		int y1 = y0 + move_y(gen);
		//		int patch_radius = radius(gen);
		//		int steps = 1 + sqrt((x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0));
		//		for (int count = 0; count < steps; ++count)
		//		{
		//			p->Dab(x0 + (x1 - x0) * count / steps, y0 + (y1 - y0) * count / steps, patch_radius, 0.8, 0.01, p_num);
		//		}
		//		std::cout << "Patch. ";
		//		try
		//		{
		//			p->Process();
		//		}
		//		catch (std::runtime_error e)
		//		{
		//			std::cout << "Exception: " << e.what() << "\n";
		//			exit(1);
		//		}
		//		catch (...)
		//		{
		//			std::cout << "Unhandled exception.\n";
		//			exit(1);
		//		}
		//	}
		//}

		
		if (false == p->Render(data))
		{
			throw std::runtime_error("Error rendering paper.\n");
		}
		

		std::stringstream ss_wcolor;
		ss_wcolor << WCOLORTEST;
		s_paper = ss_wcolor.str();

		if (0 == stbi_write_png(s_paper.c_str(), width, height, 3, data, width * 3))
		{
			throw std::runtime_error("Unable to write out wcolor_test image.\n");
		}


		//delete(pgmnt);
		delete(p);
		free(data);
		//delete(brush);
	}
	catch (std::runtime_error e)
	{
		std::cout << "Exception: Main() " << e.what() << "\n";
		exit(1);
	}
	catch (...)
	{
		std::cout << "Unhandled exception: Main()\n";
		exit(1);
	}
}

