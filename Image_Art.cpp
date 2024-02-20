// Copyright (c) 2023-2024 Steve Young
// Licensed under the MIT License

#include "Image_Art.h"
//#include <chrono>

#define PREWINDOW 3
#define POSTWINDOW 0
#define POSTGAP 0
#define XDIV 100 
#define YDIV 100
#define TEST
#define REPEAT 2
#define K 0
#define COLORSIMILARITY 17
#define PRE_ERODE_SHAPE 1
#define PRE_DILATE_SHAPE 1
#define POST_ERODE_SHAPE 1
#define POST_DILATE_SHAPE 1
#define CHANNEL 0
#define NCHANNEL 0
#define DIAGONALS 0
#define GRADTHICKNESS 3
#define PATH "D:\\temp\\"
#define OUTPUT "output"
#define FINE false;
#define POLYGON true;
#define SHOW_GRAYS false;
#define SHOW_EDGES false;
#define PALETTE 0
#define EARLY_PALETTE false;
#define SEEDS_OUT false;
#define CLOSE_FIRST false;
#define TEST_FILE "SNC00015.jpg"
// File outputs: 1: gray+edge+skeleton+paintpath, 2: base PNG, 4: base SVG, 8: post PNG, 16: post SVG, 32: paint
#define FILE_OUTPUT 255


int main(int argc, char** argv)
{
	WorkSpace* workspace = NULL;
	int repeat = REPEAT;
	int xdiv = XDIV;
	int ydiv = YDIV;
	float k = K;
	int preproc_windowsize = PREWINDOW;
	int colorsimilarity = COLORSIMILARITY;
	int pre_erode_shape = PRE_ERODE_SHAPE;
	int pre_dilate_shape = PRE_DILATE_SHAPE;
	int post_erode_shape = POST_ERODE_SHAPE;
	int post_dilate_shape = POST_DILATE_SHAPE;
	int channel = CHANNEL;
	int nchannel = NCHANNEL;
	int postproc_windowsize = POSTWINDOW;
	int postgap = POSTGAP;
	bool diagonals = DIAGONALS;
	int gradthickness = GRADTHICKNESS;
	bool fine = FINE;
	bool polygon = POLYGON;
	bool show_grays = SHOW_GRAYS;
	bool show_edges = SHOW_EDGES;
	int palette = PALETTE;
	bool early_palette = EARLY_PALETTE;
	bool seeds_out = SEEDS_OUT;
	bool close_first = CLOSE_FIRST;
	int contrast_radius = CONTRAST_RADIUS;
	unsigned char file_output = FILE_OUTPUT;

	ImageData* data_revised = NULL;

	std::string path = "";
	std::string inpath = "";
	std::string spfile = "";
	std::string grayfile = "";
	std::string edgefile = "";
	std::string temppath = "";
	std::string seeds_in = "";

	Paint_Properties prop;

	try {
		std::string filename = "";
#ifdef TEST
		filename = TEST_FILE;
		path = PATH;
#endif // TEST
		std::string tag;
		std::string value;

		int argument = 1;
		while (argument < argc)
		{
			argument++;
			std::string input = argv[argument - 1];
			int equal_loc = input.find("=");
			if (equal_loc != std::string::npos)
			{
				tag = input.substr(0, equal_loc);
				value = input.substr(equal_loc + 1, input.length() - equal_loc - 1);
			}
			else {
				tag = input;
				argument++;
				if (argument >= argc)
				{
					throw std::runtime_error("Failed to parse arguments.\n");
				}
				input = argv[argument - 1];
				if (input.find("=") == 0)
				{
					if (input.length() > 1)
					{
						value = input.substr(1, input.length() - 1);
					}
					else {
						argument++;
						if (argument >= argc)
						{
							throw std::runtime_error("Failed to parse arguments.\n");
						}
						value = argv[argument - 1];
					}
				}
				else {
					value = input;
				}
			}

			if ((tag == "b") || (tag == "background"))
			{
				prop.background = atoi(value.c_str());
			}
			else if ((tag == "bristles") || (tag == "num_bristles"))
			{
				prop.bristles = atoi(value.c_str());
			}
			else if ((tag == "brush_width") || (tag == "bw"))
			{
				prop.brush_width = atoi(value.c_str());
			}
			else if ((tag == "c") || (tag == "channel"))
			{
				channel = atoi(value.c_str());
			}
			else if ((tag == "close_first") || (tag == "reverse"))
			{
				if (atoi(value.c_str()) > 0)
				{
					close_first = true;
				}
				else {
					close_first = false;
				}
			}
			else if (tag == "colormatch")
			{
				colorsimilarity = atoi(value.c_str());
			}
			else if (tag == "contrast")
			{
				contrast_radius = atoi(value.c_str());
			}
			else if ((tag == "d") || (tag == "diagonals") || (tag == "diagonal"))
			{
				if (atoi(value.c_str()) > 0)
				{
					diagonals = true;
				}
				else {
					diagonals = false;
				}
			}
			else if ((tag == "dilate") || (tag == "predilate") || (tag == "pre_dilate") || (tag == "pre-dilate"))
			{
				pre_dilate_shape = atoi(value.c_str());
			}
			else if (tag == "edgefile")
			{
				edgefile = value.c_str();
			}
			else if ((tag == "ep") || (tag == "early_palette"))
			{
				if (atoi(value.c_str()) > 0)
				{
					early_palette = true;
				}
				else {
					early_palette = false;
				}
			}
			else if ((tag == "erode") || (tag == "preerode") || (tag == "pre_erode") || (tag == "pre-erode"))
			{
				pre_erode_shape = atoi(value.c_str());
			}
			else if ((tag == "f") || (tag == "filename"))
			{
				filename = value;
			}
			else if ((tag == "fo") || (tag == "file_output") || (tag == "output"))
			{
				file_output = atoi(value.c_str());
			}
			else if (tag == "fine")
			{
				if (atoi(value.c_str()) > 0)
				{
					fine = true;
				}
				else {
					fine = false;
				}
			}
			else if (tag == "flow")
			{
				prop.flow = atof(value.c_str());
			}
			else if ((tag == "fv") || (tag == "flow_variation"))
			{
				prop.flow_variation = atof(value.c_str());
			}
			else if ((tag == "g") || (tag == "gap"))
			{
				postgap = atof(value.c_str());
			}
			else if (tag == "glitch")
			{
				unsigned char c = atoi(value.c_str());
				if (1 == (c & 1))
				{
					prop.glitch1 = true;
				}
				else {
					prop.glitch1 = false;
				}
				if (2 == (c & 2))
				{
					prop.glitch2 = true;
				}
				else {
					prop.glitch2 = false;
				}
				if (4 == (c & 4))
				{
					prop.glitch3 = true;
				}
				else {
					prop.glitch3 = false;
				}
			}
			else if ((tag == "glitch1"))
			{
				if (atoi(value.c_str()) > 0)
				{
					prop.glitch1 = true;
				}
				else {
					prop.glitch1 = false;
				}
			}
			else if ((tag == "glitch2"))
			{
				if (atoi(value.c_str()) > 0)
				{
					prop.glitch2 = true;
				}
				else {
					prop.glitch2 = false;
				}
			}
			else if ((tag == "glitch3"))
			{
				if (atoi(value.c_str()) > 0)
				{
					prop.glitch3 = true;
				}
				else {
					prop.glitch3 = false;
				}
			}
			else if ((tag == "grad") || (tag == "gradthickness") || (tag == "grad_thickness"))
			{
				gradthickness = atoi(value.c_str());
			}
			else if (tag == "grayfile")
			{
				grayfile = value.c_str();
			}
			else if (tag == "inpath")
			{
				inpath = value.c_str();
				if (inpath.substr(inpath.length() - 1, 1) != "\\")
				{
					inpath.append("\\");
				}
			}
			else if ((tag == "k") || (tag == "snap"))
			{
				k = atof(value.c_str());
			}
			else if ((tag == "mix") || (tag == "mix_paints"))
			{
				if (atoi(value.c_str()) > 0)
				{
					prop.mix_paints = true;
				}
				else {
					prop.mix_paints = false;
				}
			}
			else if ((tag == "nc") || (tag == "negchannel") || (tag == "nchannel"))
			{
				nchannel = atoi(value.c_str());
			}
			else if ((tag == "o") || (tag == "outlines") || (tag == "outline"))
			{
				if (atoi(value.c_str()) > 0)
				{
					prop.outline = true;
				}
				else {
					prop.outline = false;
				}
			}
			else if ((tag == "p") || (tag == "post_window"))
			{
				postproc_windowsize = atoi(value.c_str());
			}
			else if ((tag == "pal") || (tag == "palette"))
			{
				palette = atoi(value.c_str());
			}
			else if (tag == "path")
			{
				path = value.c_str();
				if (path.substr(path.length() - 1, 1) != "\\")
				{
					path.append("\\");
				}
			}
			else if ((tag == "poly") || (tag == "polygon"))
			{
				if (atoi(value.c_str()) > 0)
				{
					polygon = true;
				}
				else {
					polygon = false;
				}
			}
			else if ((tag == "postdilate") || (tag == "post_dilate") || (tag == "post-dilate"))
			{
				post_dilate_shape = atoi(value.c_str());
			}
			else if ((tag == "posterode") || (tag == "post_erode") || (tag == "post-erode"))
			{
				post_erode_shape = atoi(value.c_str());
			}
			else if ((tag == "r") || (tag == "refine") || (tag == "repeat"))
			{
				repeat = atoi(value.c_str());
			}
			else if ((tag == "rv") || (tag == "radius_variation"))
			{
				if (atoi(value.c_str()) > 0)
				{
					prop.radius_variation = true;
				}
				else {
					prop.radius_variation = false;
				}
			}
			else if ((tag == "scale") || (tag == "paint_scale"))
			{
				prop.paint_scale = atof(value.c_str());
			}
			else if ((tag == "show_edges") || (tag == "edges"))
			{
				if (atoi(value.c_str()) > 0)
				{
					show_edges = true;
				}
				else {
					show_edges = false;
				}
			}
			else if ((tag == "show_grays") || (tag == "grays"))
			{
				if (atoi(value.c_str()) > 0)
				{
					show_grays = true;
				}
				else {
					show_grays = false;
				}
			}
			else if ((tag == "si") || (tag == "seeds_in"))
			{
				seeds_in = value.c_str();
			}
			else if ((tag == "so") || (tag == "seeds_out"))
			{
				if (atoi(value.c_str()) > 0)
				{
					seeds_out = true;
				}
				else {
					seeds_out = false;
				}
			}
			else if ((tag == "sp") || (tag == "subpixel") || (tag == "sub_pixel"))
			{
				if (atoi(value.c_str()) > 0)
				{
					prop.sub_pixel = true;
				}
				else {
					prop.sub_pixel = false;
				}
			}
			else if (tag == "spfile")
			{
				spfile = value.c_str();
			}
			else if ((tag == "thin") || (tag == "bristle_thin") || (tag == "thinness"))
			{
				prop.bristle_thin_factor = atof(value.c_str());
			}
			else if ((tag == "w") || (tag == "pre_window") || (tag == "windowsize"))
			{
				preproc_windowsize = atoi(value.c_str());
			}
			else if ((tag == "width_override") || (tag == "wo"))
			{
				if (atoi(value.c_str()) > 0)
				{
					prop.brush_width_override = true;
				}
				else {
					prop.brush_width_override = false;
				}
			}
			else if ((tag == "x") || (tag == "xdiv"))
			{
				xdiv = atoi(value.c_str());
			}
			else if ((tag == "y") || (tag == "ydiv"))
			{
				ydiv = atoi(value.c_str());
			}
			else {
				std::cout << "Failed to match tag value: " << tag << ", ignoring.\n";
			}
		}

		if (filename == "")
		{
			throw std::runtime_error("Need to specify filename.\n");
		}
		if (gradthickness < 3)
		{
			gradthickness = 3;
		}
		if ((postproc_windowsize - postgap) < 3)
		{
			postproc_windowsize = 0;
			postgap = 0;
		}
		std::cout << "Filename: " << filename << "\n";
		std::cout << "Show Grayscales: " << show_grays << " Show Edges: " << show_edges << "\n";
		std::cout << "SPFile: " << spfile << " Gray file: " << grayfile << " Edge file: " << edgefile << "\n";
		std::cout << "channel: " << channel << " nchannel: " << nchannel << "\n";
		std::cout << "Gradient Thickness: " << gradthickness << "\n";
		std::cout << "xdiv: " << xdiv << ", ydiv: " << ydiv << "\n";
		std::cout << "Refine: " << repeat << "\n";
		std::cout << "preprocesing window size: " << preproc_windowsize << "\n";
		std::cout << "Preprocessing erode mode: " << pre_erode_shape << ", preprocessing dilation mode: " << pre_dilate_shape << "\n";
		std::cout << "Close operation first: " << close_first << "\n";
		std::cout << "Postprocessing window size: " << postproc_windowsize << ", gap: " << postgap << "\n";
		std::cout << "Postprocessing erode mode: " << post_erode_shape << ", postprocessing dilation mode: " << post_dilate_shape << "\n";
		std::cout << "Box factor: " << k << "\n";
		std::cout << "Colormatch: " << colorsimilarity << "\n";
		std::cout << "Diagonals: " << diagonals << " Outlines: " << prop.outline << "\n";
		std::cout << "Path: " << path << " Inpath: " << inpath << "\n";
		std::cout << "Fine points for SVG: " << fine << " Polygons for SVG: " << polygon << "\n";
		std::cout << "Use contrast: " << contrast_radius << "\n";
		std::cout << "Palette size: " << palette << " Early palette reduction: " << early_palette << " Mix paints: " << prop.mix_paints << "\n";
		std::cout << "Paint scale: " << prop.paint_scale << "\n";
		std::cout << "Bristle coefficient: " << prop.bristles << " Sub-pixel: " << prop.sub_pixel << " Bristle thinness factor: " << prop.bristle_thin_factor << "\n";
		std::cout << "Brush width: " << prop.brush_width << " Radius variation: " << prop.radius_variation << " Brush width override: " << prop.brush_width_override << "\n";
		std::cout << "Flow: " << prop.flow << " Flow variation: " << prop.flow_variation << "\n";
		std::cout << "Background color: ";

		if (0 == prop.background)
		{
			std::cout << "white ";
		}
		else {
			std::cout << "black ";
		}
		std::cout << "Glitch1: ";
		if (false == prop.glitch1)
		{
			std::cout << "false";
		}
		else {
			std::cout << "true";
		}
		std::cout << " Glitch2: ";
		if (false == prop.glitch2)
		{
			std::cout << "false";
		}
		else {
			std::cout << "true";
		}
		std::cout << " Glitch3: ";
		if (false == prop.glitch3)
		{
			std::cout << "false\n";
		}
		else {
			std::cout << "true\n";
		}
		std::cout << "Seeds out: " << seeds_out << " Seeds in: " << seeds_in << "\n";
		if (show_grays || show_edges)
		{
			file_output = 1;
		}
		std::cout << "File outputs: 1: gray+edge+skeleton+paintpath, 2: base PNG, 4: base SVG, 8: post PNG, 16: post SVG, 32: paint\n";
		std::cout << "Outputs:\n";
		if (1 == (file_output & 1))
		{
			std::cout << "gray, edge, skeleton, paint path\n";
		}
		if (2 == (file_output & 2))
		{
			std::cout << "base PNG\n";
		}
		if (4 == (file_output & 4))
		{
			std::cout << "base SVG\n";
		}
		if (8 == (file_output & 8))
		{
			std::cout << "post PNG\n";
		}
		if (16 == (file_output & 16))
		{
			std::cout << "post SVG\n";
		}
		if (32 == (file_output & 32))
		{
			std::cout << "paint\n";
		}

		if ("" != inpath)
		{
			if ("" == spfile)
			{
				spfile = inpath;
				spfile.append("SuperPixels.dat");
			}
			if ("" == grayfile)
			{
				grayfile = inpath;
				grayfile.append(OUTPUT).append("_gray.png");
			}
			if ("" == edgefile)
			{
				edgefile = inpath;
				edgefile.append(OUTPUT).append("_edge.png");
			}
		}

		if (show_grays)
		{
			unsigned char* data = NULL;
			int width, height, colorchannels;
			ImageData* image = NULL;
			GradData* gray = NULL;
			data = stbi_load(filename.c_str(), &width, &height, &colorchannels, 0);
			if (NULL != data)
			{
				std::cout << "X: " << width << " Y: " << height << " Colorchannels: " << colorchannels << "\n";
			}
			else {
				throw (std::runtime_error("Unable to load image.\n"));
			}

			image = new ImageData(data, width, height, colorchannels);
			if (NULL == image)
			{
				throw (std::runtime_error("Failed to create ImageData object.\n"));
			}
			std::cout << "show_grays mode.\n";
			gray = image->gen_gray(0, 0);
			ImageData* gray_image = gray->Gradient2Image(1);
			temppath = path;
			temppath.append(OUTPUT).append("_gray");
			gray_image->write_file(temppath);
			delete gray_image;
			delete gray;
			std::cout << ".";
			gray = image->gen_gray(1, 0);
			gray_image = gray->Gradient2Image(1);
			temppath = path;
			temppath.append(OUTPUT).append("_red");
			gray_image->write_file(temppath);
			delete gray_image;
			delete gray;
			std::cout << ".";
			gray = image->gen_gray(2, 0);
			gray_image = gray->Gradient2Image(1);
			temppath = path;
			temppath.append(OUTPUT).append("_green");
			gray_image->write_file(temppath);
			delete gray_image;
			delete gray;
			std::cout << ".";
			gray = image->gen_gray(3, 0);
			gray_image = gray->Gradient2Image(1);
			temppath = path;
			temppath.append(OUTPUT).append("_blue");
			gray_image->write_file(temppath);
			delete gray_image;
			delete gray;
			std::cout << ".";
			gray = image->gen_gray(1, 2);
			gray_image = gray->Gradient2Image(1);
			temppath = path;
			temppath.append(OUTPUT).append("_red_green");
			gray_image->write_file(temppath);
			delete gray_image;
			delete gray;
			std::cout << ".";
			gray = image->gen_gray(1, 3);
			gray_image = gray->Gradient2Image(1);
			temppath = path;
			temppath.append(OUTPUT).append("_red_blue");
			gray_image->write_file(temppath);
			delete gray_image;
			delete gray;
			std::cout << ".";
			gray = image->gen_gray(2, 3);
			gray_image = gray->Gradient2Image(1);
			temppath = path;
			temppath.append(OUTPUT).append("_green_blue");
			gray_image->write_file(temppath);
			delete gray_image;
			delete gray;
			std::cout << ".\n";
			delete image;
			std::cout << "Done.\n";
		}
		else if (show_edges)
		{
			unsigned char* data = NULL;
			int width, height, colorchannels;
			ImageData* image = NULL;
			GradData* gray = NULL;
			GradData* preprocessed_gray = NULL;
			GradData* edge;

			data = stbi_load(filename.c_str(), &width, &height, &colorchannels, 0);
			if (NULL != data)
			{
				std::cout << "X: " << width << " Y: " << height << " Colorchannels: " << colorchannels << "\n";
			}
			else {
				throw (std::runtime_error("Unable to load image.\n"));
			}

			image = new ImageData(data, width, height, colorchannels);
			if (NULL == image)
			{
				throw (std::runtime_error("Failed to create ImageData object.\n"));
			}
			std::cout << "show_edges mode.\n";
			gray = image->gen_gray(channel, nchannel);
			ImageData* gray_image = gray->Gradient2Image(1);
			temppath = path;
			temppath.append(OUTPUT).append("_gray");
			gray_image->write_file(temppath);
			delete gray_image;
			std::cout << ".";
			int steps;
			int modes = 0;

			if (close_first)
			{
				steps = 6;
			}
			else {
				steps = 9;
			}
			if (pre_erode_shape)
			{
				modes += steps;
			}
			if (pre_dilate_shape)
			{
				modes += (15 - steps);
			}

			int w_size;
			if (preproc_windowsize < 13)
			{
				preproc_windowsize = 13;
			}
			for (w_size = 3; w_size <= preproc_windowsize; w_size += (int)(((preproc_windowsize - 3) / 5)))
			{
				preprocessed_gray = gray->Preprocess_Gray(4, steps, modes, w_size);
				edge = preprocessed_gray->Generate_Gradient(1, gradthickness, xdiv, ydiv, k);
				gray_image = edge->Gradient2Image(0);
				temppath = path;
				temppath.append(OUTPUT).append("_edge_").append(std::to_string(w_size));
				gray_image->write_file(temppath);
				delete gray_image;
				delete edge;
				delete preprocessed_gray;
				std::cout << ".";
			}

			delete gray;
			delete image;
			std::cout << "\nDone.\n";
		}
		else if ("" != spfile)
		{
			std::cout << "Reading SuperPixel input from file.\n";
			workspace = new WorkSpace(spfile, diagonals);
			if (NULL == workspace)
			{
				throw std::runtime_error("Failed to create Workspace.\n");
			}
			std::cout << ".";
			if (early_palette)
			{
				if (false == workspace->ReduceToPalette(0, palette))
				{
					throw std::runtime_error("Error reducing color palette.\n");
				}
				palette = 0; // Avoid re-creating palette later.
				std::cout << ".";
			}
			workspace->FindPaths(0, polygon, fine);
			if (4 == (file_output & 4))
			{
				std::cout << ".";
				temppath = path;
				temppath.append("SuperPixels.svg");
				workspace->WriteSuperPixelsSVG(temppath, 0, polygon, fine, palette, contrast_radius);
				std::cout << ".";
			}

			if (colorsimilarity > 0)
			{
				for (int ri = 0; ri < 3; ri++)
				{
					std::cout << "\nStarting Absorb run " << (ri + 1) << " ";
					if (false == workspace->CombineSuperPixels(colorsimilarity))
					{
						throw std::runtime_error("Error combining SuperPixels.\n");
					}
				}
			}
			if (postproc_windowsize > 0)
			{
				std::cout << "\nStart postprocessing ";
				workspace->Postprocess_SuperPixels(1, 1, post_erode_shape, postproc_windowsize);
				std::cout << ".";
				workspace->Postprocess_SuperPixels(1, 0, post_dilate_shape, (postproc_windowsize - postgap));
				std::cout << ".";
				workspace->FindPaths(1, polygon, fine);
				std::cout << ".";
			}
			std::cout << "\nWriting output files ";
			if (colorsimilarity > 0)
			{
				workspace->FindPaths(0, polygon, fine);
				std::cout << ".";
				if (4 == (file_output & 4))
				{
					temppath = path;
					temppath.append("SuperPixels_b.svg");
					workspace->WriteSuperPixelsSVG(temppath, 0, polygon, fine, palette, contrast_radius);
					std::cout << ".";
				}
			}
			if (2 == (file_output & 2))
			{
				data_revised = workspace->GenerateImage(0, prop);
				temppath = path;
				temppath.append(OUTPUT);
				data_revised->write_file(temppath);
				std::cout << ".";
				delete data_revised;
			}
			if ((1 == (file_output & 1)) || (32 == (file_output & 32)))
			{
				workspace->ThinSuperPixels(prop.glitch3);
				if (1 == (file_output & 1))
				{
					data_revised = workspace->GenerateImage(2, prop);
					temppath = path;
					temppath.append(OUTPUT).append("_skeleton");
					data_revised->write_file(temppath);
					std::cout << ".";
					delete data_revised;
				}
			}
			if (32 == (file_output & 32))
			{
				data_revised = workspace->GenerateImage(3, prop);
				temppath = path;
				temppath.append(OUTPUT).append("_paint");
				data_revised->write_file(temppath);
				std::cout << ".";
				delete data_revised;
			}
			if (1 == (file_output & 1))
			{
				temppath = path;
				temppath.append("SuperPixels_Paint_Paths.svg");
				workspace->WritePaintCurvesSVG(temppath);
				std::cout << ".";
			}
			if (postproc_windowsize > 0)
			{
				if (16 == (file_output & 16))
				{
					temppath = path;
					temppath.append("SuperPixels_Post.svg");
					workspace->WriteSuperPixelsSVG(temppath, 1, polygon, fine, palette, contrast_radius);
					std::cout << ".";
				}
				if (8 == (file_output & 8))
				{
					data_revised = workspace->GenerateImage(1, prop);
					temppath = path;
					temppath.append(OUTPUT).append("_post");
					data_revised->write_file(temppath);
					std::cout << ".";
					delete data_revised;
				}
			}
			std::cout << "\nDone.\n";
		}
		else {
			std::cout << "Starting from original image.";
			std::cout << "\nGenerating grayscale and gradient\n";
			workspace = new WorkSpace(filename, channel, nchannel, diagonals);
			if (NULL == workspace)
			{
				throw std::runtime_error("Failed to create Workspace.\n");
			}
			std::cout << ".";
			if (1 == (file_output & 1))
			{
				ImageData* gray_image = workspace->Gradient2Image(4);
				temppath = path;
				temppath.append(OUTPUT).append("_gray");
				gray_image->write_file(temppath);
				delete gray_image;
				std::cout << ".";
			}
			int steps;
			int modes = 0;

			if (close_first)
			{
				steps = 6;
			}
			else {
				steps = 9;
			}
			if (pre_erode_shape)
			{
				modes += steps;
			}
			if (pre_dilate_shape)
			{
				modes += (15 - steps);
			}

			workspace->Preprocess_Gray(4, steps, modes, preproc_windowsize);
			std::cout << ".";
			workspace->Generate_Gradient(1, gradthickness, xdiv, ydiv, k);
			std::cout << ".";
			if (1 == (file_output & 1))
			{
				ImageData* gradient_image = workspace->Gradient2Image(0);
				temppath = path;
				temppath.append(OUTPUT).append("_edge");
				gradient_image->write_file(temppath);
				delete gradient_image;
				std::cout << ".";
			}

			if ("" != seeds_in)
			{
				temppath = path;
				temppath.append(seeds_in);
				workspace->InitialSuperPixels(temppath);
			}
			else {
				workspace->InitialSuperPixels();
			}

			if (seeds_out)
			{
				temppath = path;
				temppath.append("seeds.txt");
				workspace->WriteSeeds(temppath);
				std::cout << ".";
			}

			while (repeat > 0)
			{
				std::cout << "\nBeginning watershed calculation.  Repeat count: " << repeat;
				repeat--;
				if (false == workspace->Watershed())
				{
					throw std::runtime_error("Error calling Watershed function.\n");
				}
				if (false == workspace->SetAveColors())
				{
					throw std::runtime_error("Error setting average colors for the SuperPixels.\n");
				}

				if (0 == repeat)
				{
					if (path.size() > 0)
					{
						temppath = path;
						temppath.append("SuperPixels.dat");
						workspace->WriteSuperPixels(temppath);
						std::cout << ".";
					}
					if (colorsimilarity > 0)
					{
						for (int ri = 0; ri < 3; ri++)
						{
							std::cout << "\nStarting Absorb run " << (ri + 1) << " ";
							if (false == workspace->CombineSuperPixels(colorsimilarity))
							{
								throw std::runtime_error("Error combining SuperPixels.\n");
							}
							if (false == workspace->SetAveColors())
							{
								throw std::runtime_error("Error setting average colors for the SuperPixels.\n");
							}
						}
					}
					else {
						if (false == workspace->SetAveColors())
						{
							throw std::runtime_error("Error setting average colors for the SuperPixels.\n");
						}
					}
					if (postproc_windowsize > 0)
					{
						std::cout << "\nStart postprocessing ";
						workspace->Postprocess_SuperPixels(1, 1, post_erode_shape, postproc_windowsize);
						std::cout << ".";
						workspace->Postprocess_SuperPixels(1, 0, post_dilate_shape, (postproc_windowsize - postgap));
						std::cout << ".";
						workspace->FindPaths(1, polygon, fine);
						std::cout << ".";
					}
				}

				if (2 == (file_output & 2))
				{
					// Paint each pixel in each SuperPixel to the average color.
					data_revised = workspace->GenerateImage(0, prop);
					temppath = path;
					temppath.append(OUTPUT);
					data_revised->write_file(temppath);
					delete data_revised;
					std::cout << ".";
				}
				if (0 == repeat)
				{
					std::cout << "\nWriting output files ";
					if (seeds_out)
					{
						temppath = path;
						temppath.append("seeds.txt");
						workspace->WriteSeeds(temppath);
						std::cout << ".";
					}
					if ((1 == (file_output & 1)) || (32 == (file_output & 32)))
					{
						workspace->ThinSuperPixels(prop.glitch3);
						if (1 == (file_output & 1))
						{
							data_revised = workspace->GenerateImage(2, prop);
							temppath = path;
							temppath.append(OUTPUT).append("_skeleton");
							data_revised->write_file(temppath);
							delete data_revised;
							std::cout << ".";
						}
					}
					if (early_palette)
					{
						if (false == workspace->ReduceToPalette(0, palette))
						{
							throw std::runtime_error("Error reducing color palette.\n");
						}
						palette = 0; // Avoid re-creating palette.
						std::cout << ".";
					}
					workspace->FindPaths(0, polygon, fine);
					if (32 == (file_output & 32))
					{
						data_revised = workspace->GenerateImage(3, prop);
						temppath = path;
						temppath.append(OUTPUT).append("_paint");
						data_revised->write_file(temppath);
						delete data_revised;
						std::cout << ".";
					}
					if (1 == (file_output & 1))
					{
						temppath = path;
						temppath.append("SuperPixels_Paint_Paths.svg");
						workspace->WritePaintCurvesSVG(temppath);
						std::cout << ".";
					}
					if (postproc_windowsize > 0)
					{
						if (16 == (file_output & 16))
						{
							temppath = path;
							temppath.append("SuperPixels_Post.svg");
							workspace->WriteSuperPixelsSVG(temppath, 1, polygon, fine, palette, contrast_radius);
							std::cout << ".";
						}
						if (8 == (file_output & 8))
						{
							data_revised = workspace->GenerateImage(1, prop);
							temppath = path;
							temppath.append(OUTPUT).append("_post");
							data_revised->write_file(temppath);
							std::cout << ".";
							delete data_revised;
						}
					}
					if (4 == (file_output & 4))
					{
						workspace->FindPaths(0, polygon, fine);
						std::cout << ".";
						temppath = path;
						temppath.append("SuperPixels.svg");
						workspace->WriteSuperPixelsSVG(temppath, 0, polygon, fine, palette, contrast_radius);
						std::cout << ".";
					}
				}

				if (repeat > 0)
				{
					if (false == workspace->SplitSuperPixels(2.0))
					{
						throw std::runtime_error("Error splitting SuperPixels.\n");
					}
				}

			}
			std::cout << "\nDone.\n";
		}
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
}




