#include "Paper.h"
#include "Brush_CUDA.cuh"

bool Paper::PaperSlope()  // Induces velocities in the water due to paper slopes.
{
	bool ret = true;
	ret = SlopeVelocities(thickness, u, v, w, h, gauss_radius);
	if (!ret)
	{
		throw std::runtime_error("Call to SlopeVelocities failed.\n");
	}
	return ret;
}

bool Paper::Calc_M_prime(int x, int y)
{
	bool ret = true;
	int x0, x1, y0, y1;
	if ((x < 0) || (y < 0))
	{
		x0 = 0;
		y0 = 0;
		x1 = w - 1;
		y1 = h - 1;
	}
	else {
		x0 = std::max(0, x - K_radius);
		y0 = std::max(0, y - K_radius);
		x1 = std::min(w - 1, x + K_radius);
		y1 = std::min(h - 1, y + K_radius);
	}
	if (!CalcMprime(M, M_prime, M_prime_kernel, K_radius, M_prime_kernel_total, w, h, x0, y0, x1, y1))
	{
		throw std::runtime_error("Failed in Calc_M_prime.\n");
		ret = false;
	}
	return ret;
}

bool Paper::OutputWaterConcentration(int layer)
{
	bool ret = true;
	std::string name;
	std::stringstream ss;
	name.clear();
	name.append("g");
	if (layer >= 0)
	{
		ss << std::setw(4) << std::setfill('0') << layer;
		name.append(ss.str());
	}
	name.append(".png");
	ret = WriteOutSparseFloatMatrix(pigments[layer]->Get_g(), w, h, name, 0.0f, 5.0f);
	return ret;
}

bool Paper::OutputDeposition(int layer)
{
	bool ret = true;
	std::string name;
	std::stringstream ss;
	name.clear();
	name.append("d");
	if (layer >= 0)
	{
		ss << std::setw(4) << std::setfill('0') << layer;
		name.append(ss.str());
	}
	name.append(".png");
	ret = WriteOutSparseFloatMatrix(pigments[layer]->Get_d(), w, h, name, 0.0f, 5.0f);
	return ret;
}

bool Paper::OutputPressure(int frame)
{
	bool ret = true;
	std::string name;
	std::stringstream ss;
	name.clear();
	name.append("p");
	if (frame > 0)
	{
		ss << std::setw(2) << std::setfill('0') << frame;
		name.append(ss.str());
	}
	name.append(".png");
	ret = WriteOutFloatArray(p, w, h, name, -0.25f, 5.0f);
	return ret;
}

bool Paper::OutputSaturation(int frame)
{
	bool ret = true;
	std::string name;
	std::stringstream ss;
	name.clear();
	name.append("s");
	if (frame > 0)
	{
		ss << std::setw(2) << std::setfill('0') << frame;
		name.append(ss.str());
	}
	name.append(".png");
	ret = WriteOutFloatArray(s, w, h, name, 0.0f, 1.0f);
	return ret;
}

bool Paper::Render(unsigned char* data)
{
	bool ret = true;
	float substrate_color[3];
	for (int channel_index = 0; channel_index < 3; ++channel_index)
	{
		substrate_color[channel_index] = color.channel[channel_index]/255.0f;
	}
	ret = CRender(pigments, data, w, h, substrate_color);
	if (!ret)
	{
		throw std::runtime_error("Error in Render.\n");
	}
	return ret;
}

bool Paper::Dry()
{
	bool ret = true;
	int num_pigments = pigments.size();
	SparseFloatMatrix* local_d = NULL;
	SparseFloatMatrix* local_g = NULL;

	int num_chunks = x_chunks * y_chunks;
	for (int k = 0; k < num_pigments; ++k)
	{
		dried.insert(k);
		ret = ret && CDry(pigments[k]->Get_g(), pigments[k]->Get_d());
	}
	if (!ret)
	{
		throw std::runtime_error("Error in CDry function.\n");
	}

	ret = ret && ResetFloatArray(u, w + 1, h, 0.0f);
	ret = ret && ResetFloatArray(v, w, h + 1, 0.0f);
	ret = ret && ResetFloatArray(p, w, h, 0.0f);
	ret = ret && ResetFloatArray(s, w, h, saturation_dry_value);
	ret = ret && ResetBoolArray(M, w, h, false);

	ret = ret && Calc_M_prime();
	return ret;
}

bool Paper::CheckDry(int layer)
{
	bool ret = (dried.find(layer) != dried.end());
	return ret;
}


bool Paper::SetVelocity(int x, int y, float u_vel, float v_vel, bool sum, bool wait)
{
	bool ret = true;
	ret = CSetVelocity(u, v, x, y, w, h, u_vel, v_vel, sum, wait);
	if (!ret)
	{
		throw std::runtime_error("Failed to set velocity.\n");
	}
	return ret;
}

bool Paper::SetBackgroundColor(Color c)
{
	bool ret = true;
	color = c;
	return ret;
}

float Paper::MaxVelocity()  // Find the absolute value of the maximum velocity.
{
	float ret = 0.0f;
	ret = GetMaxVelocity(u, v, w, h, results1, results2, results_length);
	return ret;
}

Paper::Paper(int width, int height, Color c, float saturation)
{
	color = c;
	float* temp_thickness = NULL;
	float gauss_thin_factor = 0.5; // 0.5
	float* gauss_kernel = NULL;
	float gauss_total = 0.0;
	float paper_average = 0.5; // .7
	float paper_range = 0.5;
	int num_fibers = width * height / 3600;
	float fiber_length = 50.0;
	int num_blobs = width * height / 2400;
	float blob_small = 3.0;
	float blob_large = 28.0;
	float* local_thickness = NULL;

	w = width;
	h = height;
	int N = w * h;
	int padded_N = (w + 2 * gauss_radius) * (h + 2 * gauss_radius);
	pigments.clear();
	dried.clear();
	x_chunks = 0;
	y_chunks = 0;
	results_length = (N + threads_per_block - 1) / threads_per_block;

	// Confirm that chunk_size is a multiple of AVX2_stride.
	if (0 != (chunk_size * chunk_size) % threads_per_block)
	{
		throw std::runtime_error("Value of chunk_size is not a multiple of block size.\n");
		return;
	}

	if (threads_per_block != block_dimension * block_dimension)
	{
		throw std::runtime_error("threads_per_block is not the square of block_dimension.\n");
		return;
	}

	// Paper array initializations
	M = BoolArray(w, h, false);
	thickness = FloatArray(w + 2 * gauss_radius, h + 2 * gauss_radius, true, paper_average);
	g = FloatArray(w, h, true);
	g_prime = FloatArray(w, h, true);
	u = FloatArray(w + 1, h, true);
	v = FloatArray(w, h + 1, true);
	p = FloatArray(w, h, true);
	s = FloatArray(w, h, true, saturation);
	s_prime = FloatArray(w, h, true);
	results1 = FloatArray(results_length, 1, false);
	results2 = FloatArray(results_length, 1, false);
	u_prime = FloatArray(w + 1, h, true);
	v_prime = FloatArray(w, h + 1, true);
	M_prime = FloatArray(w, h, true);
	M_prime_kernel = GaussKernel(K_radius, gauss_thin_factor);
	M_prime_kernel_total = GaussTotal(M_prime_kernel, K_radius);
	delta_matrix = FloatArray(w, h, false);

	x_chunks = (w + chunk_size - 1) / chunk_size;
	y_chunks = (h + chunk_size - 1) / chunk_size;
	M_chunks = BoolArray(x_chunks, y_chunks, false);
	int num_chunks = x_chunks * y_chunks;
	host_M_chunks = (bool*)malloc(sizeof(bool) * num_chunks);
	if (NULL == host_M_chunks)
	{
		throw std::runtime_error("Failed to allocate memory for host_M_chunks.\n");
		return;
	}

	// temp_thickness initialization
	temp_thickness = FloatArray(w + 2 * gauss_radius, h + 2 * gauss_radius, false);

	local_thickness = (float*)malloc(padded_N * sizeof(float));
	if (NULL == local_thickness)
	{
		throw std::runtime_error("Failed to allocate memory for local_thickness in Paper::Paper.\n");
		return;
	}

	std::random_device rd;
	std::mt19937 gen(rd());

	std::uniform_real_distribution<> thick_dist(paper_average - paper_range, paper_average + paper_range);

	for (int pos = 0; pos < padded_N; ++pos)
	{
		local_thickness[pos] = thick_dist(gen);
	}

	// Initialization of gauss_kernel and gauss_total.
	gauss_kernel = GaussKernel(gauss_radius, gauss_thin_factor);
	gauss_total = GaussTotal(gauss_kernel, gauss_radius);

	// Use the gauss_kernel to smooth the random thickness of the paper.
	if (!CopyFromHost(local_thickness, padded_N, temp_thickness))
	{
		throw std::runtime_error("Failed to copy memory for local_thickness to device in Paper::Paper.\n");
		return;
	}

	if (!Convolve(temp_thickness, w, h, gauss_radius, thickness, gauss_kernel, gauss_total, gauss_radius, paper_average))
	{
		throw std::runtime_error("Call to Convolve failed in Paper::Paper.\n");
		return;
	}

	if (!CopyToHost(thickness, padded_N, local_thickness))
	{
		throw std::runtime_error("Failed to copy memory for local_thickness from device in Paper::Paper.\n");
		return;
	}

	// New gauss_kernel and gauss_total.
	FreeGaussKernel(gauss_kernel);
	gauss_thin_factor = gauss_thin_factor * 0.5;
	gauss_kernel = GaussKernel(gauss_radius, gauss_thin_factor);
	gauss_total = GaussTotal(gauss_kernel, gauss_radius);

	//Add random fiber lines.
	std::uniform_real_distribution<> fb_length(fiber_length * 0.8, fiber_length * 1.2);
	std::uniform_real_distribution<> fb_rotation(0.0f, M_PI / 4.0);
	std::uniform_int_distribution<> x_rnd(0, width - 1);
	std::uniform_int_distribution<> y_rnd(0, height - 1);
	paper_range = paper_range / 2.5;
	for (int i = 0; i < num_fibers; ++i)
	{
		int x = x_rnd(gen);
		int y = y_rnd(gen);
		float fiber_value = (float)thick_dist(gen);
		float length = fb_length(gen);
		if (0 == i % 2)  // Horizontal
		{
			float dy = sin(fb_rotation(gen));
			if (0 == i % 4)
			{
				dy = -dy;
			}
			float wx = length * sqrt(1.0 - dy * dy);
			float wx_end = x - wx;
			wx = wx + x;
			float wy = length * dy + y;
			while (wx >= wx_end)
			{
				if ((wx >= 0) && (wx < width) && (wy >= 0) && (wy < height))
				{
					local_thickness[gauss_radius + (int)wx + ((int)wy + gauss_radius) * (w + 2 * gauss_radius)] = fiber_value;
				}
				wx -= 1.0;
				wy -= dy;
			}
		}
		else {  // Vertical
			float dx = sin(fb_rotation(gen));
			if (0 == (i + 1) % 4)
			{
				dx = -dx;
			}
			float wy = length * sqrt(1.0 - dx * dx);
			float wy_end = y - wy;
			wy = wy + y;
			float wx = length * dx + x;
			while (wy >= wy_end)
			{
				if ((wx >= 0) && (wx < width) && (wy >= 0) && (wy < height))
				{
					local_thickness[gauss_radius + (int)wx + ((int)wy + gauss_radius) * (w + 2 * gauss_radius)] = fiber_value;
				}
				wy -= 1.0;
				wx -= dx;
			}
		}
	}

	// Add blobs
	std::uniform_real_distribution<> blob_radius(blob_small, blob_large);
	for (int i = 0; i < num_blobs; ++i)
	{
		int x = x_rnd(gen);
		int y = y_rnd(gen);
		float radius = blob_radius(gen);
		float blob_value = 0.0625 * ((float)thick_dist(gen) - paper_average);
		for (int wx = x - radius; wx <= x + radius; ++wx)
		{
			for (int wy = y - radius; wy <= y + radius; ++wy)
			{
				if ((wx >= 0) && (wx < width) && (wy >= 0) && (wy < height))
				{
					float dx = x - wx;
					float dy = y - wy;
					if ((dx * dx + dy * dy) <= (radius * radius))
					{
						local_thickness[gauss_radius + (int)wx + ((int)wy + gauss_radius) * (w + 2 * gauss_radius)] += blob_value;
					}
				}
			}
		}
	}

	// Use the gauss_kernel to smooth the random thickness of the paper.
	if (!CopyFromHost(local_thickness, padded_N, temp_thickness))
	{
		throw std::runtime_error("Failed to copy memory for local_thickness to device in Paper::Paper.\n");
		return;
	}

	if (!Convolve(temp_thickness, w, h, gauss_radius, thickness, gauss_kernel, gauss_total, gauss_radius, paper_average))
	{
		throw std::runtime_error("Call to Convolve failed in Paper::Paper.\n");
		return;
	}

	if (!RenormalizeFloatArray(thickness, w + 2 * gauss_radius, h + 2 * gauss_radius, paper_average - paper_range, paper_average + paper_range))
	{
		throw std::runtime_error("Failed to renormalize thickness in Paper constructor.\n");
	}

	// Free up memory allocated for use only in this function.
	FreeGaussKernel(gauss_kernel);
	FreeFloatArray(temp_thickness);
	free(local_thickness);
}

Paper::~Paper()
{
	FreeGaussKernel(M_prime_kernel);
	FreeFloatArray(M_prime);
	FreeFloatArray(g);
	FreeFloatArray(g_prime);
	FreeFloatArray(v_prime);
	FreeFloatArray(u_prime);
	FreeFloatArray(s_prime);
	FreeFloatArray(results1);
	FreeFloatArray(results2);
	FreeFloatArray(s);
	FreeFloatArray(p);
	FreeFloatArray(v);
	FreeFloatArray(u);
	FreeFloatArray(thickness);
	FreeBoolArray(M);
	FreeFloatArray(delta_matrix);
	FreeBoolArray(M_chunks);
	if (NULL != host_M_chunks)
	{
		free(host_M_chunks);
	}

	for (std::vector<Pigment*>::iterator pig_it = pigments.begin(); pig_it != pigments.end(); ++pig_it)
	{
		Pigment* pigment_pntr = *pig_it;
		if (NULL != pigment_pntr)
		{
			delete(pigment_pntr);
		}
	}
	pigments.clear();
}

int Paper::GetWidth()
{
	return w;
}

int Paper::GetHeight()
{
	return h;
}

bool* Paper::GetM()
{
	return M;
}

float* Paper::GetS()
{
	return s;
}

float* Paper::GetP()
{
	return p;
}

float* Paper::GetFullG()
{
	return g;
}

int Paper::GetPadding()
{
	return gauss_radius;
}

bool Paper::GetThickness(float* out_matrix)
{
	// Take the thickness data from the GPU and put it (without padding) into the out_matrix.
	bool ret = true;
	if (NULL == out_matrix)
	{
		throw std::runtime_error("Paper::GetThickness called with NULL pointer.\n");
		ret = false;
	}
	if (ret)
	{
		ret = PartialCopyToHost(thickness, w + 2 * gauss_radius, h + 2 * gauss_radius, gauss_radius, gauss_radius, out_matrix, w, h, 0, 0, w, h);
	}
	return ret;
}

bool Paper::Tick(bool slope_velocity, bool debug_output)
{
	// Set slope_velocity to true in order for paper slopes to induce velocity changes.
	bool ret = true;
	try
	{
		if (debug_output)
		{
			OutputPressure(1);
		}
		ret = MoveWater(slope_velocity, debug_output);
		if (!ret)
		{
			throw std::runtime_error("MoveWater() failed.\n");
		}
		if (debug_output)
		{
			OutputPressure(6);
		}
		if (ret)
		{
			ret = MovePigment();
			if (!ret)
			{
				throw std::runtime_error("MovePigment() failed.\n");
			}
		}
		if (debug_output)
		{
			OutputPressure(7);
		}
		if (ret)
		{
			ret = TransferPigment();
			if (!ret)
			{
				throw std::runtime_error("TransferPigment() failed.\n");
			}
		}
		if (ret)
		{
			ret = CapillaryFlow();
			if (!ret)
			{
				throw std::runtime_error("CapillaryFlow() failed.\n");
			}
		}
		if (debug_output)
		{
			OutputPressure(8);
		}
	}
	catch (std::runtime_error e)
	{
		std::cout << "Exception: Tick() " << e.what() << "\n";
		exit(1);
	}
	catch (...)
	{
		std::cout << "Unhandled exception: Tick()\n";
		exit(1);
	}
	return ret;
}

std::vector<Pigment*> Paper::GetPigments()
{
	return pigments;
}

int Paper::SetPigment(Pigment* pgmnt)
{
	int ret = 0;
	if (NULL == pgmnt)
	{
		throw std::runtime_error("Null pointer passed to SetPigment.\n");
		ret = -1;
	}
	if (ret >= 0)
	{
		pigments.push_back(pgmnt);
		ret = pigments.size() - 1;
	}
	return ret;
}

bool Paper::Dab(int x, int y, int radius, float saturation, float concentration, int pgmnt_num)
{
	bool ret = true;
	if ((pgmnt_num >= 0) && (pgmnt_num < pigments.size()))
	{
		Pigment* pgmnt = pigments[pgmnt_num];
		ret = CalcDab(M, s, p, pgmnt->Get_g(), w, h, x, y, radius, saturation, concentration);
		if (!ret)
		{
			throw std::runtime_error("Failed to apply pigment in Dab.\n");
		}
		pgmnt->Get_d()->SyncChunks(pgmnt->Get_g());
	}
	dried.erase(pgmnt_num);
	return ret;
}

bool Paper::PaintArea(int* SPdata, int w, int h, RectQuad window, int identifier, float saturation, float concentration, int pgmnt_num)
{
	bool ret = true;
	if ((pgmnt_num >= 0) && (pgmnt_num < pigments.size()))
	{
		Pigment* pgmnt = pigments[pgmnt_num];
		ret = CPaintArea(M, s, p, pgmnt->Get_g(), SPdata, w, h, window, identifier, saturation, concentration);
		if (!ret)
		{
			throw std::runtime_error("Failed to apply pigment in Paper::PaintArea.\n");
		}
		pgmnt->Get_d()->SyncChunks(pgmnt->Get_g());
	}
	dried.erase(pgmnt_num);
	return ret;
}

Pigment::Pigment(int width, int height, float Kr_in, float Kg_in, float Kb_in, float Sr_in, float Sg_in, float Sb_in, float rho_in, float omega_in, float gamma_in)
{
	w = width;
	h = height;
	K[0] = Kr_in;
	K[1] = Kg_in;
	K[2] = Kb_in;
	S[0] = Sr_in;
	S[1] = Sg_in;
	S[2] = Sb_in;
	rho = rho_in;
	omega = omega_in;
	gamma = gamma_in;

	g = new SparseFloatMatrix(w, h, 0.0f);
	if (NULL == g)
	{
		throw std::runtime_error("Failed to create SparseFloatMatrix for g in Pigment.\n");
		return;
	}
	d = new SparseFloatMatrix(w, h, 0.0f);
	if (NULL == d)
	{
		throw std::runtime_error("Failed to create SparseFloatMatrix for d in Pigment.\n");
		return;
	}
	x_chunks = g->GetXChunks();
	y_chunks = g->GetYChunks();
}

Pigment::Pigment(int width, int height, Color* c, float r, float rho_in, float omega_in, float gamma_in)
{
	w = width;
	h = height;
	UpdateColor(c, r);
	rho = rho_in;
	omega = omega_in;
	gamma = gamma_in;

	g = new SparseFloatMatrix(w, h, 0.0f);
	if (NULL == g)
	{
		throw std::runtime_error("Failed to create SparseFloatMatrix for g in Pigment.\n");
		return;
	}
	d = new SparseFloatMatrix(w, h, 0.0f);
	if (NULL == d)
	{
		throw std::runtime_error("Failed to create SparseFloatMatrix for d in Pigment.\n");
		return;
	}
	x_chunks = g->GetXChunks();
	y_chunks = g->GetYChunks();
}

Pigment::~Pigment()
{
	if (NULL != g)
	{
		delete g;
	}
	if (NULL != d)
	{
		delete d;
	}
}

SparseFloatMatrix* Pigment::Get_g()
{
	return g;
}

SparseFloatMatrix* Pigment::Get_d()
{
	return d;
}

float Pigment::Get_rho()
{
	return rho;
}

float Pigment::Get_omega()
{
	return omega;
}

float Pigment::Get_gamma()
{
	return gamma;
}

int Pigment::GetXChunks()
{
	return x_chunks;
}

int Pigment::GetYChunks()
{
	return y_chunks;
}

bool Pigment::UpdateColor(Color* c, float r)
{
	// Set the S and K values (per Kubelka-Munk model) for the pigment to approximate the desired color.
	// c is the color being approximated when painted over a white surface.
	// r is the ratio of the reflectivity when painted over a black surface compared to a white surface.  Should be in the range of 0.1 to 0.9.

	bool ret = true;
	if (r > 0.9f)
	{
		r = 0.9f;
	}
	else if (r < 0.1f)
	{
		r = 0.1f;
	}
	for (int color_chan = 0; color_chan < 3; ++color_chan)
	{
		float Rw = (float)c->channel[color_chan] / 256.0; // Reflectivity on a white surface.  Used 256 instead of 255 to ensure the value is less than 1.0.
		if (Rw < 0.02f)
		{
			Rw = 0.02f;
		}
		float Rb = Rw * r;
		if (Rb < 0.01f)
		{
			Rb = 0.01f;
		}
		float a = 0.5 * (Rw + (Rb - Rw + 1.0f) / Rb);
		float b = sqrt(a * a - 1.0f);
		S[color_chan] = arccoth((b * b - (a - Rw) * (a - 1.0f)) / (b * (1.0f - Rw))) / b;
		K[color_chan] = S[color_chan] * (a - 1);
	}
	return ret;
}

float Pigment::GetK(int channel)
{
	float ret = 0.0f;
	if ((channel >= 0) && (channel < 3))
	{
		ret = K[channel];
	}
	return ret;
}

float Pigment::GetS(int channel)
{
	float ret = 0.0f;
	if ((channel >= 0) && (channel < 3))
	{
		ret = S[channel];
	}
	return ret;
}

bool Paper::UpdateM_chunks()
{
	bool ret = CalcActiveChunks(M, M_chunks, host_M_chunks, x_chunks, y_chunks, w, h);
	if (!ret)
	{
		throw std::runtime_error("Error in CalcActiveChunks.\n");
	}
	return ret;
}

bool Paper::Process(bool out)
{
	bool ret = UpdateM_chunks();
	if (!ret)
	{
		throw std::runtime_error("UpdateM_chunks failed in Process.\n");
	}
	else {
		if (out)
		{
			int num_pgmnt = pigments.size();
			for (int pgmnt_index = 0; pgmnt_index < num_pgmnt; ++pgmnt_index)
			{
				OutputWaterConcentration(pgmnt_index);
				OutputDeposition(pgmnt_index);
			}
			OutputPressure(1);
		}
		ret = ret && Calc_M_prime();
		std::cout << "\n";
		for (int count = 0; ret && (count <= process_steps); ++count)
		{
			try {
				ret = Tick(count < slope_velocity_steps, out);
				if (!ret)
				{
					throw std::runtime_error("Error in Tick called from Process.\n");
				}
				int num_pgmnt = pigments.size();
				std::cout << ".";
			}
			catch (std::runtime_error e)
			{
				std::cout << "Exception: " << e.what() << "\n";
				std::string explanation = e.what();
				throw std::runtime_error(explanation);
			}
			catch (...)
			{
				std::cout << "Unhandled exception in Process.\n";
				std::string explanation = "Unhandled error in Tick called from Process.\n";
				throw std::runtime_error(explanation);
			}
		}
	}

	Dry();
	std::cout << " Done.\n";
	return ret;
}

bool Paper::MoveWater(bool slope_velocity, bool debug_output)
{
	bool ret = true;

	if (debug_output)
	{
		OutputPressure(2);

		std::string name;
		name.clear();
		name.append("M_prime.png");
		ret = WriteOutFloatArray(M_prime, w, h, name, 0.0f, 1.0f);
		name.clear();
		name.append("u0.png");
		ret = WriteOutFloatArray(u, w + 1, h, name, -0.05f, 0.05f);
		name.clear();
		name.append("v0.png");
		ret = ret && WriteOutFloatArray(v, w, h + 1, name, -0.05f, 0.05f);
	}


	ret = ret && UpdateVelocities(slope_velocity);

	if (debug_output)
	{
		OutputPressure(3);

		std::string name;
		name.clear();
		name.append("u1.png");
		ret = WriteOutFloatArray(u, w + 1, h, name, -0.05f, 0.05f);
		name.clear();
		name.append("v1.png");
		ret = ret && WriteOutFloatArray(v, w, h + 1, name, -0.05f, 0.05f);
	}

	ret = ret && RelaxDivergence();

	if (debug_output)
	{
		OutputPressure(4);

		std::string name;
		name.clear();
		name.append("u2.png");
		ret = WriteOutFloatArray(u, w + 1, h, name, -0.05f, 0.05f);
		name.clear();
		name.append("v2.png");
		ret = ret && WriteOutFloatArray(v, w, h + 1, name, -0.05f, 0.05f);
	}


	ret = ret && FlowOutward();

	if (debug_output)
	{
		OutputPressure(5);

		std::string name;
		name.clear();
		name.append("u3.png");
		ret = WriteOutFloatArray(u, w + 1, h, name, -0.05f, 0.05f);
		name.clear();
		name.append("v3.png");
		ret = ret && WriteOutFloatArray(v, w, h + 1, name, -0.05f, 0.05f);
	}
	return ret;
}

bool Paper::MovePigment()  // *** Consider some mixing between cells not based on velocity. ***
{
	// Strategy:
	// In the Paper constructor, allocate two scratch variables (g, g_prime) to be used in the MoveWater calculations.  These are full matrices for the entire w x h of the
	// image, not instances of SparseMatrix.  These are zeroed out (CUDA kernel), then populated from the SparseMatrix instances of each pigment in turn.  
	// MoveWater operates on these full matrices, and keeps tabs on any writes that go over into any previously un-allocated parts of the SparseMatrix's.  
	// This will be tracked in a vector of booleans that is of length x_chunks*y_chunks.  It is initially set to reflect which chunks are active for a 
	// pigment (when g and g_prime are populated), and then used by the host-based function to expand as necessary any SparseMatrix allocations as 
	// it is writing the updated values back to the SparseMatrix.
	bool ret = true;
	float dt = 1.0f; // Intermediate step size, delta t.
	float max_vel = MaxVelocity();
	if (max_vel > 1.0f)
	{
		dt = 1.0f / max_vel;
	}
	ret = CMovePigment(w, h, u, v, M, host_M_chunks, pigments, dt, g, g_prime);
	if (!ret)
	{
		throw std::runtime_error("Error calling CMovePigment from MovePigment.\n");
	}
	return ret;
}

bool Paper::TransferPigment()
{
	bool ret = true;
	int num_chunks = x_chunks * y_chunks;
	if (pigments.size() > 0) // If no pigments, no reason to go further.
	{
		for (std::vector<Pigment*>::iterator pgmnt_it = pigments.begin(); pgmnt_it != pigments.end(); ++pgmnt_it)
		{
			Pigment* k_pntr = *pgmnt_it;
			CTransferPigment(k_pntr, thickness, gauss_radius, M, M_chunks);
		}
	}
	return ret;
}

bool Paper::CapillaryFlow()
{
	bool ret = true;

	ret = CCapillaryFlow(w, h, s, s_prime, p, thickness, gauss_radius, M, M_chunks);
	if (!ret)
	{
		throw std::runtime_error("Error in CCapillaryFlow called from CapillaryFlow.\n");
	}
	ret = ret && Calc_M_prime();
	if (!ret)
	{
		throw std::runtime_error("Error in Calc_M_prime called from CapillaryFlow.\n");
	}
	ret = ret && UpdateM_chunks();
	if (!ret)
	{
		throw std::runtime_error("Error in UpdateM_chunks called from CapillaryFlow.\n");
	}
	return ret;
}

bool Paper::UpdateVelocities(bool slope_velocity)
{
	bool ret = true;

	float dt = 1.0f; // Intermediate step size, delta t.
	float t = 0.0f; // To keep track of the progress through the full time step.
	if (slope_velocity)
	{
		ret = PaperSlope();
	}
	float max_vel = MaxVelocity();
	if (max_vel > 1.0f)
	{
		dt = 1.0f / max_vel;
	}
	float step_size = dt; // The step size, which is also the scale factor for u and v, so that they take into account the time step size.

	ret = CalcVelocities(dt, w, h, u, v, u_prime, v_prime, p, M);

	return ret;
}

bool Paper::RelaxDivergence()
{
	bool ret = true;
	ret = CalcRelaxDivergence(M, u, v, p, delta_matrix, w, h);
	if (!ret)
	{
		throw std::runtime_error("Failed in RelaxDivergence.\n");
	}
	return ret;
}

bool Paper::FlowOutward()
{
	bool ret = true;
	if (!CalcFlowOutward(M, M_prime, p, w, h))
	{
		throw std::runtime_error("FlowOutwad failed.\n");
		ret = false;
	}
	return ret;
}
