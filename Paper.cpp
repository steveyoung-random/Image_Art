#include "Paper.h"

bool Paper::EnforceBoundaries()
{
	bool ret = true;

	// First, ensure all borders have zero velocity.
	// Start with v:
	for (int i = 0; i < w; ++i)
	{
		v[i][0] = 0.0f;
		v[i][h] = 0.0f;
	}
	// Then u:
	for (int j = 0; j < h; ++j)
	{
		u[0][j] = 0.0f;
		u[w][j] = 0.0f;
	}

	for (std::set<int>::iterator chunk_it = M_chunks.begin(); chunk_it != M_chunks.end(); ++chunk_it)
	{
		int chunk_index = *chunk_it;
		int wy = (chunk_index / x_chunks);  // Y component of beginning of window represented by chunk.
		int wx = (chunk_index % x_chunks); // X component of the beginning of windows represented by chunk. 
		bool last_row = (wy == (y_chunks - 1));
		bool last_column = (wx == (x_chunks - 1));
		bool first_row = (0 == wy);
		bool first_column = (0 == wx);
		wx = wx * chunk_size;
		wy = wy * chunk_size;
		int i_limit = std::min(wx + chunk_size, w);
		int j_limit = std::min(wy + chunk_size, h);
		bool current_mask;
		// First, manage all of the u values, by row:
		for (int j = wy; j < j_limit; ++j)
		{
			// Set up run over the row at j.
			if (!first_column)
			{
				current_mask = M[wx - 1][j];
			}
			else {
				current_mask = false; // Reflects that all areas beyond the matrix are zero.
			}
			for (int i = wx; i < i_limit; ++i)
			{
				if (M[i][j] != current_mask) // Indicates there has been a transition of mask states.
				{
					u[i][j] = 0.0f;
					current_mask = M[i][j];
				}
			}
			if (!last_column)  // This has to be done because the next chunk may not be in the mask set.
			{
				if (M[i_limit][j] != current_mask)
				{
					u[i_limit][j] = 0.0f;
				}
			}
		}
		// Next, manage all of the v values, by column:
		for (int i = wx; i < i_limit; ++i)
		{
			// Set up run over the column at i.
			if (!first_row)
			{
				current_mask = M[i][wy - 1];
			}
			else {
				current_mask = false; // Reflects that all areas beyond the matrix are zero.
			}
			for (int j = wy; j < j_limit; ++j)
			{
				if (M[i][j] != current_mask) // Indicates there has been a transition of mask states.
				{
					v[i][j] = 0.0f;
					current_mask = M[i][j];
				}
			}
			if (!last_row)  // This has to be done because the next chunk may not be in the mask set.
			{
				if (M[i][j_limit] != current_mask)
				{
					v[i][j_limit] = 0.0f;
				}
			}
		}
	}
	return ret;
}

bool Paper::PaperSlope()  // Induces velocities in the water due to paper slopes.
{
	bool ret = true;
	for (std::set<int>::iterator chunk_it = M_chunks.begin(); chunk_it != M_chunks.end(); ++chunk_it)
	{
		int chunk_index = *chunk_it;
		int wy = (chunk_index / x_chunks) * chunk_size;  // Y component of beginning of window represented by chunk.
		int wx = (chunk_index % x_chunks) * chunk_size; // X component of the beginning of windows represented by chunk. 
		int i_limit = std::min(wx + chunk_size, w);
		int j_limit = std::min(wy + chunk_size, h);
		for (int i = wx; i < i_limit; ++i)
		{
			if (M_column[chunk_index][i - wx] > 0)
			{
				for (int j = wy; j < j_limit; ++j)
				{
					if (M[i][j])
					{
						if ((M[i - 1][j]) && (i > 0)) // Only update for places not on the edge of the paper.
						{
							u[i][j] -= slope_factor * (thickness[i][j] - thickness[i - 1][j]);
						}
						if ((M[i][j - 1]) && (j > 0)) // Only update for places not on the edge of the paper.
						{
							v[i][j] -= slope_factor * (thickness[i][j] - thickness[i][j - 1]);
						}
					}
				}
			}
		}
	}
	return ret;
}

bool Paper::Calc_M_prime(int x, int y)
{
	bool ret = true;
	int limit_x0, limit_x1, limit_y0, limit_y1;
	if ((x < 0) || (y < 0))
	{
		limit_x0 = 0;
		limit_y0 = 0;
		limit_x1 = w - 1;
		limit_y1 = h - 1;
	}
	else {
		limit_x0 = std::max(0, x - K_radius);
		limit_x1 = std::min(w - 1, x + K_radius);
		limit_y0 = std::max(0, y - K_radius);
		limit_y1 = std::min(h - 1, y + K_radius);
	}
	for (int i = limit_x0; i <= limit_x1; ++i)
	{
		for (int j = limit_y0; j <= limit_y1; ++j)
		{
			float mask_accumulation = 0.0f;
			for (int wx = i - K_radius; wx <= i + K_radius; ++wx) // Window for kernel centered on i,j
			{
				for (int wy = j - K_radius; wy <= j + K_radius; ++wy)
				{
					if ((wx >= 0) && (wx < w) && (wy >= 0) && (wy < h))
					{
						if (M[wx][wy])
						{
							mask_accumulation += M_prime_kernel[abs(wx - i)][abs(wy - j)];
						}
					}
				}
			}
			M_prime[i][j] = 1.0f - (mask_accumulation / M_prime_kernel_total);
		}
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
	ret = WriteOutFloatArray(p, w, h, name, 0.0f, 5.0f);
	return ret;
}

bool Paper::Render(unsigned char* data)
{
	bool ret = true;
	float R, T;
	float substrate_color[3];
	for (int channel_index = 0; channel_index < 3; ++channel_index)
	{
		substrate_color[channel_index] = (float)(color.channel[channel_index] / 255.0f);
	}
	int num_pgmnt = pigments.size();
	for (int i = 0; i < w; ++i)
	{
		for (int j = 0; j < h; ++j)
		{
			long pos = 3 * (i + j * w);
			for (int channel_index = 0; channel_index < 3; ++channel_index)
			{
				float reflected = 0.0;
				float upward_transmission = 1.0;
				float downward_transmission = 1.0;
				for (int pgmnt_index = num_pgmnt - 1; (pgmnt_index >= 0) && (upward_transmission > EFFECTIVE_ZERO); --pgmnt_index)
				{
					Pigment* pgmnt_ptr = pigments[pgmnt_index];
					float pigment_thickness = std::max(0.0f, pgmnt_ptr->Get_d()->Get_Value(i, j) + pgmnt_ptr->Get_g()->Get_Value(i, j) * render_g_value);
					R = pigments[pgmnt_index]->GetR(channel_index, pigment_thickness);
					T = pigments[pgmnt_index]->GetT(channel_index, pigment_thickness);
					reflected += downward_transmission * R * upward_transmission;
					// Set up for next layer.
					downward_transmission = downward_transmission * (1.0 - R) * T;
					upward_transmission = upward_transmission * T;
				}
				// White substrate.
				reflected += downward_transmission * upward_transmission * substrate_color[channel_index];
				data[pos + channel_index] = 255 * reflected;
			}
		}
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
		if (0.0f != render_g_value)
		{
			local_d = pigments[k]->Get_d();
			local_g = pigments[k]->Get_g();
			for (int chunk_index = 0; chunk_index < num_chunks; ++chunk_index)
			{
				if (local_g->CheckChunkNumber(chunk_index))
				{
					if (!local_d->CheckChunkNumber(chunk_index))
					{
						int wy = chunk_index / x_chunks;  // Y component of beginning of window represented by chunk.
						int wx = chunk_index - (wy * x_chunks); // X component of the beginning of windows represented by chunk. 
						wx = wx * chunk_size;
						wy = wy * chunk_size;
						local_d->Set_Value(wx, wy, 0.0f, false); // Initialize the chunk.
					}
					float** g_chunk = local_g->GetChunk(chunk_index);
					float** d_chunk = local_d->GetChunk(chunk_index);
					for (int wx = 0; wx < chunk_size; ++wx)
					{
						for (int wy = 0; wy < chunk_size; ++wy)
						{
							d_chunk[wx][wy] += render_g_value * g_chunk[wx][wy];
						}
					}
				}
			}
		}
	}

	for (int i = 0; i <= w; ++i)
	{
		for (int j = 0; j <= h; ++j)
		{
			if (j < h)
			{
				u[i][j] = 0.0;
			}
			if (i < w)
			{
				v[i][j] = 0.0;
			}
			if ((i < w) && (j < h))
			{
				M[i][j] = false;
				p[i][j] = 0.0f;
				s[i][j] = std::min(s[i][j], saturation_dry_value);
			}
		}
	}
	for (int k = 0; k < num_pigments; ++k)
	{
		ret = ret && pigments[k]->Get_g()->Reset_Values();
	}
	for (int chunk_id = 0; chunk_id < (x_chunks * y_chunks); ++chunk_id)
	{
		for (int i = 0; i < chunk_size; ++i)
		{
			M_column[chunk_id][i] = 0; // Reset all mask column values to zero.
		}
	}
	ret = ret && Calc_M_prime();
	return ret;
}

bool Paper::CheckDry(int layer)
{
	bool ret = (dried.find(layer) != dried.end());
	return ret;
}

bool Paper::SetVelocity(int x, int y, float u_vel, float v_vel, bool sum)
{
	bool ret = true;
	if ((x >= 0) && (x < w) && (y >= 0) && (y < h))
	{
		if (sum)
		{
			u[x][y] += u_vel;
			u[x + 1][y] += u_vel;
			v[x][y] += v_vel;
			v[x][y + 1] += v_vel;
		}
		else {
			u[x][y] = u_vel;
			u[x + 1][y] = u_vel;
			v[x][y] = v_vel;
			v[x][y + 1] = v_vel;
		}
	}
	return ret;
}

bool Paper::SetBackgroundColor(Color c)
{
	bool ret = true;
	color = c;
	return ret;
}

float Paper::MaxVelocity()  // Find the absolute value of the maximum velocity within all chunks with non-zero mask.
{
	float ret = 0.0f;
	float value = 0.0f;
	for (std::set<int>::iterator chunk_it = M_chunks.begin(); chunk_it != M_chunks.end(); ++chunk_it)
	{
		int chunk_index = *chunk_it;
		int wy = (chunk_index / x_chunks) * chunk_size;  // Y component of beginning of window represented by chunk.
		int wx = (chunk_index % x_chunks) * chunk_size; // X component of the beginning of windows represented by chunk. 
		bool last_row = (wy == (y_chunks - 1));
		bool last_column = (wx == (x_chunks - 1));
		int i_limit = std::min(wx + chunk_size, w);
		int j_limit = std::min(wy + chunk_size, h);
		for (int i = wx; i < i_limit; ++i)
		{
			for (int j = wy; j < j_limit; ++j)
			{
				value = abs(u[i][j]);
				if (value > ret)
				{
					ret = value;
				}
				value = abs(v[i][j]);
				if (value > ret)
				{
					ret = value;
				}
			}
			if (last_row)
			{
				value = abs(v[i][h]);
				if (value > ret)
				{
					ret = value;
				}
			}
		}
		if (last_column)
		{
			for (int j = wy; j < j_limit; ++j)
			{
				value = abs(u[w][j]);
				if (value > ret)
				{
					ret = value;
				}
			}
		}
	}
	return ret;
}

float Paper::capacity(int i, int j)
{
	float ret = thickness[i][j] * (capacity_max - capacity_min) + capacity_min;
	return ret;
}

float Paper::u_vel(int i, int j)
{
	float ret = 0.0f;
	if ((i >= 0) && (i < w) && (j >= 0) && (j < h))
	{
		ret = (u[i][j] + u[i + 1][j]) / 2.0f;
	}
	return ret;
}

float Paper::v_vel(int i, int j)
{
	float ret = 0.0f;
	if ((i >= 0) && (i < w) && (j >= 0) && (j < h))
	{
		ret = (v[i][j] + v[i][j + 1]) / 2.0f;
	}
	return ret;
}

float Paper::uv_corner(int i, int j, int x, int y)
{
	// x, y values are 0 or 1.  0 means -0.5 and 1 means +0.5 for the x or y dimension.  The combination of x and y values determine which corner is being selected.
	// i, j values select the location of the center of the cell.
	// Return value is the product of the u and v values multipled together, each as calculated at the selected corner.

	float ret = 0.0f;
	float u_val, v_val; // Sum of u and v components, which each need to be divided by 2 to get average values.
	if ((i >= 0) && (i < w) && (j >= 0) && (j < h) && (x >= 0) && (x <= 1) && (y >= 0) && (y <= 1))
	{
		if (j + y == 0)
		{
			u_val = u[i + x][0];
		}
		else if (j + y == h)
		{
			u_val = u[i + x][j - 1 + y];
		}
		else {
			u_val = (u[i + x][j - 1 + y] + u[i + x][j + y]);
		}
		if (i + x == 0)
		{
			v_val = v[0][j + y];
		}
		else if (i + x == w)
		{
			v_val = v[i - 1 + x][j + y];
		}
		else {
			v_val = (v[i - 1 + x][j + y] + v[i + x][j + y]);
		}
		ret = u_val * v_val / 4.0f;
	}
	return ret;
}

Paper::Paper(int width, int height, Color c, float saturation)
{
	color = c;
	float** temp_thickness = NULL;
	int gauss_radius = 50; // 50
	float gauss_thin_factor = 0.5; // 0.5
	float** gauss_kernel = NULL;
	float gauss_total = 0.0;
	float paper_average = 0.5; // .7
	float paper_range = 0.5;
	int num_fibers = width * height / 1200;
	float fiber_length = 20.0;
	int num_blobs = width * height / 2400;
	float blob_small = 3.0;
	float blob_large = 8.0;

	long pos; // Position in the one-dimenisional data structure.
	w = width;
	h = height;
	pigments.clear();
	dried.clear();
	x_chunks = 0;
	y_chunks = 0;
	M_column = NULL;

	// Confirm that chunk_size is a multiple of AVX2_stride.
	if (0 != chunk_size % AVX2_stride)
	{
		throw std::runtime_error("Value of chunk_size is not a multiple of AVX2_stride.\n");
		return;
	}


	// Paper array initializations
	M = BoolArray(w, h, false);
	M_chunks.clear();
	thickness = FloatArray(w, h, false);
	u = FloatArray(w + 1, h, true);
	v = FloatArray(w, h + 1, true);
	p = FloatArray(w, h, true);
	s = FloatArray(w, h, true, saturation);
	s_prime = FloatArray(w, h, true);
	g = new SparseFloatMatrix(w, h, 0.0f);
	if (NULL == g)
	{
		throw std::runtime_error("Failed to allocate SparseFloatMatrix for g in Paper constructor.\n");
		return;
	}
	g_prime = new SparseFloatMatrix(w, h, 0.0f);
	if (NULL == g_prime)
	{
		throw std::runtime_error("Failed to allocate SparseFloatMatrix for g_prime in Paper constructor.\n");
		return;
	}
	u_prime = FloatArray(w + 1, h, true);
	v_prime = FloatArray(w, h + 1, true);
	M_prime = FloatArray(w, h, true);
	M_prime_kernel = GaussKernel(K_radius, gauss_thin_factor);
	M_prime_kernel_total = GaussTotal(M_prime_kernel, K_radius);

	x_chunks = g->GetXChunks();
	y_chunks = g->GetYChunks();
	int num_chunks = x_chunks * y_chunks;
	u_delta = (float***)malloc(sizeof(float**) * num_chunks);
	if (NULL == u_delta)
	{
		throw std::runtime_error("Failed to allocate memory for u_delta.\n");
		return;
	}
	v_delta = (float***)malloc(sizeof(float**) * num_chunks);
	if (NULL == v_delta)
	{
		throw std::runtime_error("Failed to allocate memory for v_delta.\n");
		return;
	}
	M_column = (int**)malloc(sizeof(int*) * num_chunks);
	if (NULL == M_column)
	{
		throw std::runtime_error("Failed to allocate memory for M_column.\n");
		return;
	}
	for (int chunk_id = 0; chunk_id < num_chunks; ++chunk_id)
	{
		u_delta[chunk_id] = FloatArray(chunk_size + 1, chunk_size, true);
		if (NULL == u_delta[chunk_id])
		{
			throw std::runtime_error("Failed to allocate memory for u_delta[].\n");
			return;
		}
		v_delta[chunk_id] = FloatArray(chunk_size, chunk_size + 1, true);
		if (NULL == v_delta[chunk_id])
		{
			throw std::runtime_error("Failed to allocate memory for v_delta[].\n");
			return;
		}
		M_column[chunk_id] = (int*)malloc(sizeof(int) * chunk_size);
		if (NULL == M_column[chunk_id])
		{
			throw std::runtime_error("Failed to allocate memory for M_column[chunk_id].\n");
			return;
		}
		for (int i = 0; i < chunk_size; ++i)
		{
			M_column[chunk_id][i] = 0;
		}
	}



	// temp_thickness initialization
	temp_thickness = FloatArray(w, h, false);


	std::random_device rd;
	std::mt19937 gen(rd());

	std::uniform_real_distribution<> thick_dist(paper_average - paper_range, paper_average + paper_range);
	for (int i = 0; i < w; ++i)
	{
		for (int j = 0; j < h; ++j)
		{
			temp_thickness[i][j] = thick_dist(gen);
		}
	}

	// Initialization of gauss_kernel and gauss_total.
	gauss_kernel = GaussKernel(gauss_radius, gauss_thin_factor);
	gauss_total = GaussTotal(gauss_kernel, gauss_radius);

	// Use the gauss_kernel to smooth the random thickness of the paper.
	if (use_AVX2)
	{
		int padWidth, padHeight;
		float* ConInput = PadMatrix(temp_thickness, w, h, gauss_radius, padWidth, padHeight, paper_average);
		float* ConOutput = PadMatrix(temp_thickness, w, h, gauss_radius, padWidth, padHeight);
		bool ret = Convolve(ConInput, ConOutput, gauss_kernel, padWidth, padHeight, gauss_radius);
		if (!ret)
		{
			throw std::runtime_error("Convolve failed.\n");
			return;
		}
		for (int i = 0; i < w; ++i)
		{
			for (int j = 0; j < h; ++j)
			{
				thickness[i][j] = ConOutput[padHeight * (i + gauss_radius) + j + gauss_radius] / gauss_total;
			}
		}
		FreePadMatrix(ConOutput);
		FreePadMatrix(ConInput);
	}
	else {
		for (int y = 0; y < height; ++y)
		{
			for (int x = 0; x < width; ++x)
			{
				float cumulative_value = 0.0;
				for (int j = -gauss_radius; j <= gauss_radius; ++j)
				{
					for (int i = -gauss_radius; i <= gauss_radius; ++i)
					{
						int wx, wy; // x and y values iterating over the window of the kernel overlaid with the center at x,y.
						wx = x + i;
						wy = y + j;
						if ((wx < 0) || (wy < 0) || (wx >= width) || (wy >= height))
						{
							cumulative_value += paper_average * gauss_kernel[abs(i)][abs(j)];
						}
						else {
							cumulative_value += temp_thickness[wx][wy] * gauss_kernel[abs(i)][abs(j)];
						}
					}
				}
				thickness[x][y] = cumulative_value / gauss_total;
			}
		}
	}

	CopyFloatArray(thickness, temp_thickness, w, h);

	// New gauss_kernel and gauss_total.
	FreeGaussKernel(gauss_kernel, gauss_radius);
	gauss_thin_factor = gauss_thin_factor * 0.75;
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
					//temp_thickness[(int)wx][(int)wy] = std::min(temp_thickness[(int)wx][(int)wy], paper_average - paper_range);
					temp_thickness[(int)wx][(int)wy] = fiber_value;
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
					//temp_thickness[(int)wx][(int)wy] = std::min(temp_thickness[(int)wx][(int)wy], paper_average - paper_range);
					temp_thickness[(int)wx][(int)wy] = fiber_value;
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
		for (int wx = x - radius; wx <= x + radius; ++wx)
		{
			for (int wy = y - radius; wy <= y + radius; ++wy)
			{
				if ((wx >= 0) && (wx < width) && (wy >= 0) && (wy < height))
				{
					float dx = x - wx;
					float dy = y - wy;
					float blob_value = (float)thick_dist(gen);
					if ((dx * dx + dy * dy) <= (radius * radius))
					{
						//temp_thickness[(int)wx][(int)wy] = std::min(paper_average, std::min(temp_thickness[(int)wx][(int)wy], blob_value));
						temp_thickness[(int)wx][(int)wy] = blob_value;
					}
				}
			}
		}
	}

	// Use the gauss_kernel to smooth the random thickness of the paper.
	if (use_AVX2)
	{
		int padWidth, padHeight;
		float* ConInput = PadMatrix(temp_thickness, w, h, gauss_radius, padWidth, padHeight, paper_average);
		float* ConOutput = PadMatrix(temp_thickness, w, h, gauss_radius, padWidth, padHeight);
		bool ret = Convolve(ConInput, ConOutput, gauss_kernel, padWidth, padHeight, gauss_radius);
		if (!ret)
		{
			throw std::runtime_error("Convolve failed.\n");
			return;
		}
		for (int i = 0; i < w; ++i)
		{
			for (int j = 0; j < h; ++j)
			{
				thickness[i][j] = ConOutput[padHeight * (i + gauss_radius) + j + gauss_radius] / gauss_total;
			}
		}
		FreePadMatrix(ConOutput);
		FreePadMatrix(ConInput);
	}
	else {
		for (int y = 0; y < height; ++y)
		{
			for (int x = 0; x < width; ++x)
			{
				float cumulative_value = 0.0;
				for (int j = -gauss_radius; j <= gauss_radius; ++j)
				{
					for (int i = -gauss_radius; i <= gauss_radius; ++i)
					{
						int wx, wy; // x and y values iterating over the window of the kernel overlaid with the center at x,y.
						wx = x + i;
						wy = y + j;
						if ((wx < 0) || (wy < 0) || (wx >= width) || (wy >= height))
						{
							cumulative_value += paper_average * gauss_kernel[abs(i)][abs(j)];
						}
						else {
							cumulative_value += temp_thickness[wx][wy] * gauss_kernel[abs(i)][abs(j)];
						}
					}
				}
				thickness[x][y] = cumulative_value / gauss_total;
			}
		}
	}

	if (!RenormalizeFloatArray(thickness, width, height, paper_average - paper_range, paper_average + paper_range))
	{
		throw std::runtime_error("Failed to renormalize thickness in Paper constructor.\n");
	}
	// Free up memory allocated for use only in this function.
	FreeGaussKernel(gauss_kernel, gauss_radius);
	FreeFloatArray(temp_thickness, w);
}

Paper::~Paper()
{
	FreeGaussKernel(M_prime_kernel, K_radius);
	FreeFloatArray(M_prime, w);
	FreeFloatArray(v_prime, w);
	FreeFloatArray(u_prime, w + 1);
	if (NULL != g_prime)
	{
		delete g_prime;
	}
	if (NULL != g)
	{
		delete g;
	}
	FreeFloatArray(s_prime, w);
	FreeFloatArray(s, w);
	FreeFloatArray(p, w);
	FreeFloatArray(v, w);
	FreeFloatArray(u, w + 1);
	FreeFloatArray(thickness, w);
	FreeBoolArray(M, w);
	if (NULL != M_column)
	{
		for (int i = 0; i < (x_chunks * y_chunks); ++i)
		{
			free(M_column[i]);
		}
		free(M_column);
	}
	if (NULL != u_delta)
	{
		for (int i = 0; i < (x_chunks * y_chunks); ++i)
		{
			FreeFloatArray(u_delta[i], chunk_size + 1);
		}
		free(u_delta);
	}
	if (NULL != v_delta)
	{
		for (int i = 0; i < (x_chunks * y_chunks); ++i)
		{
			FreeFloatArray(v_delta[i], chunk_size);
		}
		free(v_delta);
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

float** Paper::GetThickness()
{
	return thickness;
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
		//if (0 == x_chunks)
		//{
			//x_chunks = pgmnt->GetXChunks();
			//y_chunks = pgmnt->GetYChunks();
			//int num_chunks = x_chunks * y_chunks;
			//M_column = (int**)malloc(sizeof(int*) * num_chunks);
			//if (NULL == M_column)
			//{
			//	throw std::runtime_error("Failed to allocate memory for M_column.\n");
			//	return false;
			//}
			//for (int chunk_id = 0; chunk_id < num_chunks; ++chunk_id)
			//{
			//	M_column[chunk_id] = (int*)malloc(sizeof(int) * chunk_size);
			//	if (NULL == M_column[chunk_id])
			//	{
			//		throw std::runtime_error("Failed to allocate memory for M_column[chunk_id].\n");
			//		return false;
			//	}
			//	for (int i = 0; i < chunk_size; ++i)
			//	{
			//		M_column[chunk_id][i] = 0;
			//	}
			//}
		//}
	}
	return ret;
}

bool Paper::Dab(int x, int y, int radius, float saturation, float concentration, int pgmnt_num)
{
	bool ret = true;
	if ((pgmnt_num >= 0) && (pgmnt_num < pigments.size()))
	{
		Pigment* pgmnt = pigments[pgmnt_num];
		if (radius < 1)
		{
			if ((x > 0) && (x < (w - 1)) && (y > 0) && (y < (h - 1)))
			{
				s[x][y] += saturation;
				pgmnt->Add_g(x, y, concentration);
				M[x][y] = true;
				p[x][y] += dab_pressure;
			}
		}else{
			for (int i = x - radius; i <= (x + radius); ++i)
			{
				if ((i > 0) && (i < (w - 1)))
				{
					int offset = sqrt(radius * radius - abs(x - i) * abs(x - i));
					for (int j = y - offset; j <= (y + offset); ++j)
					{
						if ((j > 0) && (j < (h - 1)))
						{
							s[i][j] += saturation;
							pgmnt->Add_g(i, j, concentration);
							M[i][j] = true;  // Note that M_column values are not updated here.  They are updated in the Process function.
							p[i][j] += dab_pressure;
						}
					}
				}
			}
		}
	}
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

bool Pigment::Add_g(int x, int y, float concentration)
{
	bool ret = true;
	if ((x >= 0) && (x < w) && (y >= 0) && (y < h))
	{
		ret = ret && g->Set_Value(x, y, concentration, true);
	}
	return ret;
}

float Pigment::GetR(int channel, float thickness)
{
	float ret = 0.0f;
	float a;
	float b;
	float c;
	float d;
	if ((channel >= 0) && (channel < 3))
	{
		a = 1.0 + K[channel] / S[channel];
		b = sqrt(a * a - 1);
		d = b * S[channel] * thickness;
		c = a * sinh(d) + b * cosh(d);
		if (c > EFFECTIVE_ZERO)
		{
			ret = sinh(d) / c;
		}
		else {
			ret = 0.0;
		}
	}
	if (ret > 1.0)
	{
		ret = 1.0;
	}
	else if (ret < 0) {
		ret = 0.0;
	}
	return ret;
}

float Pigment::GetT(int channel, float thickness)
{
	float ret = 0.0f;
	float a;
	float b;
	float c;
	float d;
	if ((channel >= 0) && (channel < 3))
	{
		a = 1.0 + K[channel] / S[channel];
		b = sqrt(a * a - 1);
		d = b * S[channel] * thickness;
		c = a * sinh(d) + b * cosh(d);
		if (c > EFFECTIVE_ZERO)
		{
			ret = b / c;
		}
		else {
			ret = 0.0;
		}
	}
	if (ret > 1.0)
	{
		ret = 1.0;
	}
	else if (ret < 0) {
		ret = 0.0;
	}
	return ret;
}

bool Pigment::UpdateColor(Color* c, float r)
{
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

bool Paper::TestMask(int chunk_index)
{
	bool ret = false;
	int wy = chunk_index / x_chunks;  // Y component of beginning of window represented by chunk.
	int wx = chunk_index - (wy * x_chunks); // X component of the beginning of windows represented by chunk. 
	wx = wx * chunk_size;
	wy = wy * chunk_size;
	for (int i = wx; (i < (wx + chunk_size)) && (i < w); ++i)
	{
		int hits = 0; // Number of true values in this column in this chunk.
		for (int j = wy; (j < (wy + chunk_size)) && (j < h); ++j)
		{
			if (M[i][j])
			{
				ret = true;
				++hits;
			}
		}
		if (0 == hits)
		{
			M_column[chunk_index][i - wx] = 0;
		}
		else if (hits < chunk_size)
		{
			M_column[chunk_index][i - wx] = 1;
		}
		else {
			M_column[chunk_index][i - wx] = 2;
		}
	}
	return ret;
}

bool Paper::UpdateM_chunks()
{
	bool ret = true;
	M_chunks.clear();
	if (x_chunks > 0)
	{
		for (int chunk_index = 0; chunk_index < x_chunks * y_chunks; ++chunk_index)  // Update M_chunks before calculations.
		{
			if (TestMask(chunk_index))
			{
				M_chunks.insert(chunk_index);
			}
		}
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
				ret = Tick(count < slope_velocity_steps, false);
				if (!ret)
				{
					throw std::runtime_error("Error in Tick called from Process.\n");
				}
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
		name.append("u0.png");
		ret = WriteOutFloatArray(u, w + 1, h, name, -0.05f, 0.05f);
		name.clear();
		name.append("v0.png");
		ret = ret && WriteOutFloatArray(v, w, h + 1, name, -0.05f, 0.05f);
	}

	ret = ret && UpdateVelocitiesThreads(slope_velocity);

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

	if (ret)
	{
		ret = RelaxDivergenceThreads();

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
	}
	if (ret)
	{
		ret = FlowOutward();

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
	}
	return ret;
}

bool Paper::MovePigment()  // *** Consider some mixing between cells not based on velocity. ***
{
	bool ret = true;
	float dt = 1.0f; // Intermediate step size, delta t.
	float t = 0.0f; // To keep track of the progress through the full time step.
	float max_vel = MaxVelocity();
	if (max_vel > 1.0f)
	{
		dt = 1.0f / max_vel;
	}
	std::vector<Pigment*>::iterator k;
	for (k = pigments.begin(); (ret && (k != pigments.end())); ++k)
	{
		Pigment* k_pntr = *k;
		float step_size = dt; // The step size, which is also the scale factor for u and v, so that they take into account the time step size.
		float keep_going = true;
		while (ret && keep_going)  // t is updated at the end of the step.
		{
			if ((false == g->Copy(k_pntr->Get_g())) || (false == g_prime->Copy(k_pntr->Get_g())))
			{
				throw std::runtime_error("Failed to copy g values in MovePigment.\n");
				ret = false;
			}
			int x_chunks = g->GetXChunks();
			int y_chunks = g->GetYChunks();
			std::set<int> chunk_set = g->GetChunkSet();
			for (std::set<int>::iterator chunk_it = chunk_set.begin(); chunk_it != chunk_set.end(); ++chunk_it)
			{
				int chunk_index = *chunk_it;
				int wy = chunk_index / x_chunks;  // Y component of beginning of window represented by chunk.
				int wx = chunk_index - (wy * x_chunks); // X component of the beginning of windows represented by chunk. 
				wx = wx * chunk_size;
				wy = wy * chunk_size;
				for (int i = wx; (i < wx + chunk_size) && (i < w); ++i)
				{
					for (int j = wy; (j < wy + chunk_size) && (j < h); ++j)
					{
						float g_move_left = 0.0f;
						float g_move_right = 0.0f;
						float g_move_top = 0.0f;
						float g_move_bottom = 0.0f;

						if (i > 0) // Left pixel available
						{
							g_move_left = std::max(0.0f, -u[i][j]);
						}
						if (i < (w - 1)) // Right pixel available
						{
							g_move_right = std::max(0.0f, u[i + 1][j]);
						}
						if (j > 0) // Top pixel available
						{
							g_move_top = std::max(0.0f, -v[i][j]);
						}
						if (j < (h - 1)) // Bottom pixel available
						{
							g_move_bottom = std::max(0.0f, v[i][j + 1]);
						}
						float outflow_factor = g_move_left + g_move_right + g_move_top + g_move_bottom;  // Detect if total outflow would be higher than 1.0.
						if (outflow_factor <= 1.0f)
						{
							outflow_factor = g->Get_Value(i, j) * step_size;
						}
						else {
							outflow_factor = g->Get_Value(i, j) * step_size / outflow_factor;
						}
						if (outflow_factor > EFFECTIVE_ZERO)
						{
							if (i > 0)
							{
								g_prime->Set_Value(i - 1, j, g_move_left * outflow_factor, true);
							}
							if (i < (w - 1))
							{
								g_prime->Set_Value(i + 1, j, g_move_right * outflow_factor, true);
							}
							if (j > 0)
							{
								g_prime->Set_Value(i, j - 1, g_move_top * outflow_factor, true);
							}
							if (j < (h - 1))
							{
								g_prime->Set_Value(i, j + 1, g_move_bottom * outflow_factor, true);
							}
							g_prime->Set_Value(i, j, -(g_move_left + g_move_right + g_move_top + g_move_bottom) * outflow_factor, true);
						}
					}
				}
			}
			if (false == k_pntr->Get_g()->Copy(g_prime))
			{
				throw std::runtime_error("Failed to copy g values in MovePigment.\n");
				ret = false;
			}

			if (t >= 1.0f)
			{
				keep_going = false;
			}
			else {
				t += dt;
				if (t > 1.0f)  // Adjust if the step would go past t = 1.0.
				{
					step_size = dt - (t - 1.0f);
					if (step_size < EFFECTIVE_ZERO)
					{
						keep_going = false;
					}
				}
			}
		}
	}
	return ret;
}

bool Paper::TransferPigment()
{
	bool ret = true;
	int num_chunks = x_chunks * y_chunks;
	if (pigments.size() > 0) // If no pigments, no reason to go further.  This allows us to assume there is at least one pigment for the TestMask call.
	{
		for (std::set<int>::iterator chunk_it = M_chunks.begin(); chunk_it != M_chunks.end(); ++chunk_it)
		{
			int chunk_index = *chunk_it;
			int wy = (chunk_index / x_chunks) * chunk_size;  // Y component of beginning of window represented by chunk.
			int wx = (chunk_index % x_chunks) * chunk_size; // X component of the beginning of windows represented by chunk. 
			int i_limit = std::min(chunk_size, w - wx);
			int j_limit = std::min(chunk_size, h - wy);
			int k;
			for (k = 0; (ret && (k < pigments.size())); ++k)
			{
				Pigment* k_pntr = pigments[k];
				bool dry_layer = CheckDry(k);
				SparseFloatMatrix* pigment_g = k_pntr->Get_g();
				SparseFloatMatrix* pigment_d = k_pntr->Get_d();
				bool g_chunk_used = pigment_g->CheckChunkNumber(chunk_index);
				bool d_chunk_used = pigment_d->CheckChunkNumber(chunk_index);
				if (g_chunk_used || d_chunk_used)  // No need to do anything if no g or d values in this chunk.
				{
					float gamma = k_pntr->Get_gamma();
					float rho = k_pntr->Get_rho();
					float omega = k_pntr->Get_omega();
					if (!g_chunk_used)
					{
						pigment_g->Set_Value(wx, wy, 0.0f, false);  // Cause initialization of this chunk for g.
					}
					if (!d_chunk_used)
					{
						pigment_d->Set_Value(wx, wy, 0.0f, false); // Cause initialization of this chunk for d.
					}
					float** chunk_g = pigment_g->GetChunk(chunk_index);  // For efficiency, work directly on the 2D g matrix for the chunk.
					float** chunk_d = pigment_d->GetChunk(chunk_index);  // For efficiency, work directly on the 2D d matrix for the chunk.
					for (int i = 0; i < i_limit; ++i)
					{
						int iwx = i + wx;
						for (int j = 0; j < j_limit; ++j)
						{
							int jwy = j + wy;
							if (M[iwx][jwy])
							{
								float delta_down = 0.0f;
								float delta_up = 0.0f;

								delta_down = chunk_g[i][j] * (1.0f - thickness[iwx][jwy] * gamma) * rho;
								delta_up = chunk_d[i][j] * (1.0f + (thickness[iwx][jwy] - 1.0f) * gamma) * rho / omega;
								if (chunk_d[i][j] + delta_down > 1.0f)
								{
									delta_down = std::max(0.0f, 1.0f - chunk_d[i][j]);
								}
								if (chunk_g[i][j] + delta_up > 1.0f)
								{
									delta_up = std::max(0.0f, 1.0f - chunk_g[i][j]);
								}
								chunk_d[i][j] += (delta_down - delta_up);
								chunk_g[i][j] += (delta_up - delta_down);
							}
						}
					}
				}
			}

		}
	}
	return ret;
}

bool Paper::CapillaryFlow()
{
	bool ret = true;

	for (std::set<int>::iterator chunk_it = M_chunks.begin(); chunk_it != M_chunks.end(); ++chunk_it)
	{
		int chunk_index = *chunk_it;
		int wy = (chunk_index / x_chunks) * chunk_size;  // Y component of beginning of window represented by chunk.
		int wx = (chunk_index % x_chunks) * chunk_size; // X component of the beginning of windows represented by chunk. 
		int i_limit = std::min(wx + chunk_size, w);
		int j_limit = std::min(wy + chunk_size, h);
		for (int i = wx; i < i_limit; ++i)
		{
			if (M_column[chunk_index][i - wx] > 0)
			{
				for (int j = wy; j < j_limit; ++j)
				{
					if (M[i][j])
					{
						s[i][j] += std::max(0.0f, std::min(absorption_alpha, capacity(i, j) - s[i][j]));
					}
				}
			}
		}
	}
	CopyFloatArray(s, s_prime, w, h);
	for (int i = 0; i < w; ++i)
	{
		for (int j = 0; j < h; ++j)
		{
			// k and l will give the four direct neighbors to i,j.
			int k = 0;
			int l = 0;
			for (int neighbor = 0; neighbor < 4; ++neighbor)  // Step through four adjacent neighbors.  Left, Up, Right, Down.
			{
				unsigned char bit0 = neighbor & 1;
				unsigned char bit1 = (neighbor & 2) >> 1;
				k = i + (1 - bit0) * (2 * bit1 - 1);
				l = j + bit0 * (2 * bit1 - 1);
				if ((k > 0) && (k < (w - 1)) && (l > 0) && (l < (h - 1))) // Check to see if k or l is out of bounds.
				{
					//					if ((s[i][j] > saturation_epsilon) && (s[i][j] > s[k][l]) && (s[k][l] > saturation_delta))
					if ((s[i][j] > saturation_epsilon * capacity(i, j)) && (s[i][j] > s[k][l]))
					{
						float del_s = std::max(0.0f, std::min(saturation_max_diffusion, std::min(s[i][j] - s[k][l], capacity(k, l) - s[k][l]) / 4.0f));
						s_prime[i][j] = s_prime[i][j] - del_s;
						s_prime[k][l] = s_prime[k][l] + del_s;
					}
				}
			}
		}
	}
	CopyFloatArray(s_prime, s, w, h);

	for (int chunk_index = 0; chunk_index < x_chunks * y_chunks; ++chunk_index)  // Update M_chunks before calculations.
	{
		bool need_update = false; // Flag to indicate whether a call to TestMask is needed.
		int wy = chunk_index / x_chunks;  // Y component of beginning of window represented by chunk.
		int wx = chunk_index - (wy * x_chunks); // X component of the beginning of windows represented by chunk. 
		wx = wx * chunk_size;
		wy = wy * chunk_size;
		int i_limit = std::min(wx + chunk_size, w);
		int j_limit = std::min(wy + chunk_size, h);
		for (int i = wx; i < i_limit; ++i)
		{
			for (int j = wy; j < j_limit; ++j)
			{
				if ((!M[i][j]) && ((s[i][j] > saturation_sigma) || (s[i][j] >= (capacity(i, j) - 2.0 * absorption_alpha))))
				{
					M[i][j] = true;
					p[i][j] = 0.0f;
					ret = ret && Calc_M_prime(i, j);
					need_update = true;
				}
			}
		}
		if (need_update)
		{
			if (TestMask(chunk_index))
			{
				M_chunks.insert(chunk_index);
			}
		}
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
	float keep_going = true;
	//	int num_chunks = x_chunks * y_chunks;
	while (ret && keep_going)  // t is updated at the end of the step.
	{
		for (std::set<int>::iterator M_c_it = M_chunks.begin(); M_c_it != M_chunks.end(); ++M_c_it)
		{
			int chunk_index = *M_c_it;
			int wy = chunk_index / x_chunks;  // Y component of beginning of window represented by chunk.
			int wx = chunk_index - (wy * x_chunks); // X component of the beginning of windows represented by chunk. 
			wx = wx * chunk_size;
			wy = wy * chunk_size;
			for (int i = wx; (i < wx + chunk_size) && (i < w); ++i)
			{
				for (int j = wy; (j < wy + chunk_size) && (j < h); ++j)
				{
					if (M[i][j])
					{
						float A, B;
						float u_ij = u_vel(i, j);
						float u_i1j = u_vel(i + 1, j);
						A = u_ij * u_ij - u_i1j * u_i1j + uv_corner(i, j, 1, 0) - uv_corner(i, j, 1, 1);
						B = u[i][j] - 4.0 * u[i + 1][j];
						if (i < (w - 1))
						{
							B += u[i + 2][j];
						}
						if (j < (h - 1))
						{
							B += u[i + 1][j + 1];
						}
						if (j > 0)
						{
							B += u[i + 1][j - 1];
						}
						if (i < (w - 1))
						{
							u_prime[i + 1][j] = u[i + 1][j] + step_size * (A + mu * B + p[i][j] - p[i + 1][j] - kappa * u[i + 1][j]);
						}
						else {
							u_prime[i + 1][j] = u[i + 1][j] + step_size * (A + mu * B + p[i][j] - kappa * u[i + 1][j]);  // Drop out the p terms that rely on information outside the bounds.
						}
						if (u_prime[i + 1][j] > max_velocity)
						{
							u_prime[i + 1][j] = max_velocity;
						}
						else if (u_prime[i + 1][j] < -max_velocity)
						{
							u_prime[i + 1][j] = -max_velocity;
						}

						float v_ij = v_vel(i, j);
						float v_ij1 = v_vel(i, j + 1);
						A = v_ij * v_ij - v_ij1 * v_ij1 + uv_corner(i, j, 0, 1) - uv_corner(i, j, 1, 1);
						B = v[i][j] - 4.0 * v[i][j + 1];
						if (i < (w - 1))
						{
							B += v[i + 1][j + 1];
						}
						if (i > 0)
						{
							B += v[i - 1][j + 1];
						}
						if (j < (h - 1))
						{
							B += v[i][j + 2];
						}
						if (j < (h - 1))
						{
							v_prime[i][j + 1] = v[i][j + 1] + step_size * (A + mu * B + p[i][j] - p[i][j + 1] - kappa * v[i][j + 1]);
						}
						else {
							v_prime[i][j + 1] = v[i][j + 1] + step_size * (A + mu * B + p[i][j] - kappa * v[i][j + 1]);  // Drop out the p terms that rely on information outside the bounds.
						}
						if (v_prime[i][j + 1] > max_velocity)
						{
							v_prime[i][j + 1] = max_velocity;
						}
						else if (v_prime[i][j + 1] < -max_velocity)
						{
							v_prime[i][j + 1] = -max_velocity;
						}
					}
				}
			}
		}
		ret = (CopyFloatArray(u_prime, u, w + 1, h) && CopyFloatArray(v_prime, v, w, h + 1));
		ret = (ret && EnforceBoundaries());
		if (t >= 1.0f)
		{
			keep_going = false;
		}
		else {
			t += dt;
			if (t > 1.0f)  // Adjust if the step would go past t = 1.0.
			{
				step_size = dt - (t - 1.0f);
				if (step_size < EFFECTIVE_ZERO)
				{
					keep_going = false;
				}
			}
		}
	}
	return ret;
}

bool Paper::UpdateVelocitiesThreads(bool slope_velocity)
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
	float keep_going = true;
	std::vector<bool> results(M_chunks.size());
	//	int num_chunks = x_chunks * y_chunks;
	while (ret && keep_going)  // t is updated at the end of the step.
	{
		std::transform(
			std::execution::par,
			M_chunks.begin(), M_chunks.end(),
			results.begin(),
			[this, step_size](int x) {
				return UpVelThreadworker(x, step_size);
			}
		);

		for (std::set<int>::iterator chunk_it = M_chunks.begin(); chunk_it != M_chunks.end(); ++chunk_it)
		{
			int chunk_index = *chunk_it;
			int wy = (chunk_index / x_chunks);  // Y component of beginning of window represented by chunk.
			int wx = (chunk_index % x_chunks); // X component of the beginning of windows represented by chunk. 
			bool last_row = (wy == (y_chunks - 1));
			bool last_column = (wx == (x_chunks - 1));
			wx = wx * chunk_size;
			wy = wy * chunk_size;
			int source_width = std::min(wx + chunk_size, w) - wx;
			int source_height = std::min(wy + chunk_size, h) - wy;
			if (last_row)
			{
				CopyPartialFloatArray(v_prime, wx, wy, source_width, source_height + 1, v, wx, wy);
			}
			else {
				CopyPartialFloatArray(v_prime, wx, wy, source_width, source_height, v, wx, wy);
			}
			if (last_column)
			{
				CopyPartialFloatArray(u_prime, wx, wy, source_width + 1, source_height, u, wx, wy);
			}
			else {
				CopyPartialFloatArray(u_prime, wx, wy, source_width, source_height, u, wx, wy);
			}

		}
		ret = (ret && EnforceBoundaries());
		if (t >= 1.0f)
		{
			keep_going = false;
		}
		else {
			t += dt;
			if (t > 1.0f)  // Adjust if the step would go past t = 1.0.
			{
				step_size = dt - (t - 1.0f);
				if (step_size < EFFECTIVE_ZERO)
				{
					keep_going = false;
				}
			}
		}
	}
	return ret;
}

bool Paper::UpVelThreadworker(int chunk_id, float step_size)
{
	int wy = chunk_id / x_chunks;  // Y component of beginning of window represented by chunk.
	int wx = chunk_id - (wy * x_chunks); // X component of the beginning of windows represented by chunk. 
	wx = wx * chunk_size;
	wy = wy * chunk_size;
	for (int i = wx; (i < wx + chunk_size) && (i < w); ++i)
	{
		for (int j = wy; (j < wy + chunk_size) && (j < h); ++j)
		{
			if (M[i][j])
			{
				float A, B;
				float u_ij = u_vel(i, j);
				float u_i1j = u_vel(i + 1, j);
				A = u_ij * u_ij - u_i1j * u_i1j + uv_corner(i, j, 1, 0) - uv_corner(i, j, 1, 1);
				B = u[i][j] - 4.0 * u[i + 1][j];
				if (i < (w - 1))
				{
					B += u[i + 2][j];
				}
				if (j < (h - 1))
				{
					B += u[i + 1][j + 1];
				}
				if (j > 0)
				{
					B += u[i + 1][j - 1];
				}
				if (i < (w - 1))
				{
					u_prime[i + 1][j] = u[i + 1][j] + step_size * (A + mu * B + p[i][j] - p[i + 1][j] - kappa * u[i + 1][j]);
				}
				else {
					u_prime[i + 1][j] = u[i + 1][j] + step_size * (A + mu * B + p[i][j] - kappa * u[i + 1][j]);  // Drop out the p terms that rely on information outside the bounds.
				}
				if (u_prime[i + 1][j] > max_velocity)
				{
					u_prime[i + 1][j] = max_velocity;
				}
				else if (u_prime[i + 1][j] < -max_velocity)
				{
					u_prime[i + 1][j] = -max_velocity;
				}

				float v_ij = v_vel(i, j);
				float v_ij1 = v_vel(i, j + 1);
				A = v_ij * v_ij - v_ij1 * v_ij1 + uv_corner(i, j, 0, 1) - uv_corner(i, j, 1, 1);
				B = v[i][j] - 4.0 * v[i][j + 1];
				if (i < (w - 1))
				{
					B += v[i + 1][j + 1];
				}
				if (i > 0)
				{
					B += v[i - 1][j + 1];
				}
				if (j < (h - 1))
				{
					B += v[i][j + 2];
				}
				if (j < (h - 1))
				{
					v_prime[i][j + 1] = v[i][j + 1] + step_size * (A + mu * B + p[i][j] - p[i][j + 1] - kappa * v[i][j + 1]);
				}
				else {
					v_prime[i][j + 1] = v[i][j + 1] + step_size * (A + mu * B + p[i][j] - kappa * v[i][j + 1]);  // Drop out the p terms that rely on information outside the bounds.
				}
				if (v_prime[i][j + 1] > max_velocity)
				{
					v_prime[i][j + 1] = max_velocity;
				}
				else if (v_prime[i][j + 1] < -max_velocity)
				{
					v_prime[i][j + 1] = -max_velocity;
				}
			}
		}
	}
	return true;
}

bool Paper::RelaxDivergence()
{
	bool ret = true;
	int count = 0;
	float delta_max = tau + 1.0f;
	while ((delta_max > tau) && (count < relaxation_steps) && ret)
	{
		ret = ret && CopyFloatArray(u, u_prime, w + 1, h) && CopyFloatArray(v, v_prime, w, h + 1);
		delta_max = 0.0f;
		for (std::set<int>::iterator M_c_it = M_chunks.begin(); M_c_it != M_chunks.end(); ++M_c_it)
		{
			int chunk_index = *M_c_it;
			int wy = chunk_index / x_chunks;  // Y component of beginning of window represented by chunk.
			int wx = chunk_index - (wy * x_chunks); // X component of the beginning of windows represented by chunk. 
			wx = wx * chunk_size;
			wy = wy * chunk_size;
			for (int i = wx; (i < wx + chunk_size) && (i < w); ++i)
			{
				for (int j = wy; (j < wy + chunk_size) && (j < h); ++j)
				{
					if (M[i][j])
					{
						float delta = -xi * (u[i + 1][j] - u[i][j] + v[i][j + 1] - v[i][j]);
						p[i][j] += delta; // If pressure_delta is -1, multiply delta by that here.
						// Following lines removed because this clamping does not seem necessary (and affects performance).
						//if (p[i][j] < 0.0f)
						//{
						//	p[i][j] = 0.0f;
						//}

						u_prime[i + 1][j] += delta;
						u_prime[i][j] -= delta;
						v_prime[i][j + 1] += delta;
						v_prime[i][j] -= delta;
						delta_max = std::max(delta_max, abs(delta));
					}
				}
			}
		}
		ret = ret && CopyFloatArray(u_prime, u, w + 1, h) && CopyFloatArray(v_prime, v, w, h + 1);
		count++;
	}
	return ret;
}

bool Paper::RelaxDivergenceThreads()
{
	bool ret = true;
	int count = 0;
	int last_column = x_chunks - 1;
	int last_row = y_chunks - 1;
	std::vector<bool> results(M_chunks.size());
	while (count < relaxation_steps)
	{
		std::transform(
			std::execution::par,
			M_chunks.begin(), M_chunks.end(),
			results.begin(),
			[this](int x) {
				return RelDivThreadWorker(x);
			}
		);
		for (std::set<int>::iterator chunk_it = M_chunks.begin(); chunk_it != M_chunks.end(); ++chunk_it)
		{
			int chunk_index = *chunk_it;
			int wy = chunk_index / x_chunks;  // Y component of beginning of window represented by chunk.
			int wx = chunk_index % x_chunks; // X component of the beginning of windows represented by chunk. 

			bool flag_last_column = (wx == last_column);
			bool flag_last_row = (wy == last_row);
			wx = wx * chunk_size;
			wy = wy * chunk_size;
			int source_width;
			int source_height;
			if (flag_last_column)
			{
				source_width = w - wx;
			}
			else {
				source_width = chunk_size;
			}
			if (flag_last_row)
			{
				source_height = h - wy;
			}
			else {
				source_height = chunk_size;
			}
			AddPartialFloatArray(u_delta[chunk_index], source_width + 1, source_height, u, wx, wy);
			AddPartialFloatArray(v_delta[chunk_index], source_width, source_height + 1, v, wx, wy);
		}
		count++;
	}
	return ret;
}

bool Paper::RelDivThreadWorker(int chunk_index)
// In order to avoid having this thread potentially write to a value being written to by another thread (when calculating the last column),
// this function only calculates the delta values for u and v, and they are actually added to u and v for each chunk sequentially after threads all complete.
{
	bool ret = true;
	int wy = chunk_index / x_chunks;  // Y component of beginning of window represented by chunk.
	int wx = chunk_index - (wy * x_chunks); // X component of the beginning of windows represented by chunk. 
	wx = wx * chunk_size;
	wy = wy * chunk_size;
	if (use_AVX2)
	{
		__m256 zero_vector = _mm256_set1_ps(0.0f);
		__m256 neg_vector = _mm256_set1_ps(-1.0f);
		ResetFloatArrayAVX2(u_delta[chunk_index], chunk_size + 1, chunk_size, &zero_vector);
		ResetFloatArrayAVX2(v_delta[chunk_index], chunk_size, chunk_size + 1, &zero_vector);
		for (int i = wx; (i < wx + chunk_size) && (i < w); ++i)
		{
			int ii = i - wx; // The i value for the local chunk.
			if (1 == M_column[chunk_index][i - wx])
			{
				for (int j = wy; (j <= wy + chunk_size - AVX2_stride) && (j < h); j += AVX2_stride)
				{
					int jj = j - wy; // The j value for the local chunk.
					__m256 maskValues = _mm256_castsi256_ps(_mm256_setr_epi32(
						M[i][j + 0] ? -1 : 0,
						M[i][j + 1] ? -1 : 0,
						M[i][j + 2] ? -1 : 0,
						M[i][j + 3] ? -1 : 0,
						M[i][j + 4] ? -1 : 0,
						M[i][j + 5] ? -1 : 0,
						M[i][j + 6] ? -1 : 0,
						M[i][j + 7] ? -1 : 0
					));
					__m256 xi_vector = _mm256_set1_ps(-xi);
					xi_vector = _mm256_blendv_ps(zero_vector, xi_vector, maskValues); // Use M mask values to set xi_vector to either -xi or 0, depending on the mask.
					__m256 temp_vector = _mm256_load_ps(&u[i + 1][j]);  // u[i+1][j]
					__m256 temp2_vector = _mm256_load_ps(&u[i][j]);  // u[i][j]
					__m256 delta_vector = _mm256_sub_ps(temp_vector, temp2_vector); // delta = u[i+1][j] - u[i][j]

					temp_vector = _mm256_load_ps(&v[i][j + 1]); // v[i][j+1]
					delta_vector = _mm256_add_ps(delta_vector, temp_vector); // delta = u[i+1][j] - u[i][j] + v[i][j+1]

					temp_vector = _mm256_load_ps(&v[i][j]); // v[i][j]
					delta_vector = _mm256_sub_ps(delta_vector, temp_vector); // delta = u[i+1][j] - u[i][j] + v[i][j+1] - v[i][j]

					delta_vector = _mm256_mul_ps(delta_vector, xi_vector);  // delta_vector should now be -xi * (u[i + 1][j] - u[i][j] + v[i][j + 1] - v[i][j]) where M[i][j] == true, zero otherwise.

					temp_vector = _mm256_load_ps(&p[i][j]);  // p[i][j]
					temp_vector = _mm256_add_ps(temp_vector, delta_vector);  // This should be p[i][j] += delta.
					_mm256_store_ps(&p[i][j], temp_vector);

					temp_vector = _mm256_load_ps(&u_delta[chunk_index][ii + 1][jj]); // u_delta[chunk_index][i+1][j]
					temp_vector = _mm256_add_ps(temp_vector, delta_vector);  // This should be u_delta[chunk_index][i + 1][j] += delta.
					_mm256_store_ps(&u_delta[chunk_index][ii + 1][jj], temp_vector);

					temp_vector = _mm256_load_ps(&u_delta[chunk_index][ii][jj]); // u_delta[chunk_index][i][j]
					temp_vector = _mm256_sub_ps(temp_vector, delta_vector);  // This should be u_delta[chunk_index][i][j] -= delta.
					_mm256_store_ps(&u_delta[chunk_index][ii][jj], temp_vector);

					temp_vector = _mm256_load_ps(&v_delta[chunk_index][ii][jj + 1]);  // v_delta[chunk_index][i][j + 1]
					temp_vector = _mm256_add_ps(temp_vector, delta_vector);  // This should be v_delta[chunk_index][i][j + 1] += delta.
					_mm256_store_ps(&v_delta[chunk_index][ii][jj + 1], temp_vector);

					temp_vector = _mm256_load_ps(&v_delta[chunk_index][ii][jj]); // v_delta[chunk_index][i][j]
					temp_vector = _mm256_sub_ps(temp_vector, delta_vector);  // This should be v_delta[chunk_index][i][j] -= delta.
					_mm256_store_ps(&v_delta[chunk_index][ii][jj], temp_vector);
				}
			}
			else if (2 == M_column[chunk_index][i - wx]) // No need for mask.
			{
				for (int j = wy; (j <= wy + chunk_size - AVX2_stride) && (j < h); j += AVX2_stride)
				{
					int jj = j - wy; // The j value for the local chunk.
					__m256 xi_vector = _mm256_set1_ps(-xi);
					__m256 temp_vector = _mm256_load_ps(&u[i + 1][j]);  // u[i+1][j]
					__m256 temp2_vector = _mm256_load_ps(&u[i][j]);  // u[i][j]
					__m256 delta_vector = _mm256_sub_ps(temp_vector, temp2_vector); // delta = u[i+1][j] - u[i][j]

					temp_vector = _mm256_load_ps(&v[i][j + 1]); // v[i][j+1]
					delta_vector = _mm256_add_ps(delta_vector, temp_vector); // delta = u[i+1][j] - u[i][j] + v[i][j+1]

					temp_vector = _mm256_load_ps(&v[i][j]); // v[i][j]
					delta_vector = _mm256_sub_ps(delta_vector, temp_vector); // delta = u[i+1][j] - u[i][j] + v[i][j+1] - v[i][j]

					delta_vector = _mm256_mul_ps(delta_vector, xi_vector);  // delta_vector should now be -xi * (u[i + 1][j] - u[i][j] + v[i][j + 1] - v[i][j]) where M[i][j] == true, zero otherwise.

					temp_vector = _mm256_load_ps(&p[i][j]);  // p[i][j]
					temp_vector = _mm256_add_ps(temp_vector, delta_vector);  // This should be p[i][j] += delta.
					_mm256_store_ps(&p[i][j], temp_vector);

					temp_vector = _mm256_load_ps(&u_delta[chunk_index][ii + 1][jj]); // u_delta[chunk_index][i+1][j]
					temp_vector = _mm256_add_ps(temp_vector, delta_vector);  // This should be u_delta[chunk_index][i + 1][j] += delta.
					_mm256_store_ps(&u_delta[chunk_index][ii + 1][jj], temp_vector);

					temp_vector = _mm256_load_ps(&u_delta[chunk_index][ii][jj]); // u_delta[chunk_index][i][j]
					temp_vector = _mm256_sub_ps(temp_vector, delta_vector);  // This should be u_delta[chunk_index][i][j] -= delta.
					_mm256_store_ps(&u_delta[chunk_index][ii][jj], temp_vector);

					temp_vector = _mm256_load_ps(&v_delta[chunk_index][ii][jj + 1]);  // v_delta[chunk_index][i][j + 1]
					temp_vector = _mm256_add_ps(temp_vector, delta_vector);  // This should be v_delta[chunk_index][i][j + 1] += delta.
					_mm256_store_ps(&v_delta[chunk_index][ii][jj + 1], temp_vector);

					temp_vector = _mm256_load_ps(&v_delta[chunk_index][ii][jj]); // v_delta[chunk_index][i][j]
					temp_vector = _mm256_sub_ps(temp_vector, delta_vector);  // This should be v_delta[chunk_index][i][j] -= delta.
					_mm256_store_ps(&v_delta[chunk_index][ii][jj], temp_vector);
				}
			}
		}
	}
	else {
		ret = ResetFloatArray(u_delta[chunk_index], chunk_size + 1, chunk_size, 0.0f);
		ret = ret && ResetFloatArray(v_delta[chunk_index], chunk_size, chunk_size + 1, 0.0f);
		for (int i = wx; (i < wx + chunk_size) && (i < w); ++i)
		{
			for (int j = wy; (j < wy + chunk_size) && (j < h); ++j)
			{
				if (M[i][j])
				{
					float delta = -xi * (u[i + 1][j] - u[i][j] + v[i][j + 1] - v[i][j]);
					p[i][j] += delta; // If pressure_delta is -1, multiply delta by that here.
					u_delta[chunk_index][i + 1][j] += delta;
					u_delta[chunk_index][i][j] -= delta;
					v_delta[chunk_index][i][j + 1] += delta;
					v_delta[chunk_index][i][j] -= delta;
				}
			}
		}
	}
	return ret;
}



bool Paper::FlowOutward()
{
	bool ret = true;
	for (std::set<int>::iterator chunk_it = M_chunks.begin(); chunk_it != M_chunks.end(); ++chunk_it)
	{
		int chunk_index = *chunk_it;
		int wy = (chunk_index / x_chunks) * chunk_size;  // Y component of beginning of window represented by chunk.
		int wx = (chunk_index % x_chunks) * chunk_size; // X component of the beginning of windows represented by chunk. 
		int i_limit = std::min(wx + chunk_size, w);
		int j_limit = std::min(wy + chunk_size, h);
		for (int i = wx; i < i_limit; ++i)
		{
			if (M_column[chunk_index][i - wx] > 0)
			{
				for (int j = wy; j < j_limit; ++j)
				{
					if (M[i][j])
					{
						//p[i][j] = std::max((float)(p[i][j] - eta * M_prime[i][j]), 0.0f);
						p[i][j] = (float)(p[i][j] - eta * M_prime[i][j]);
					}
				}
			}
		}
	}
	return ret;
}
