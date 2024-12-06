// Copyright (c) 2023-2024 Steve Young
// Licensed under the MIT License

#include "Brush.h"
#include "SPixelData.h"
#include <random>

Bristle::Bristle(FloatPointPair o, float fd)
{
	offset = o;
	last_loc.x = 0.0;
	last_loc.y = 0.0;
	wander.x = 0.0;
	wander.y = 0.0;
	flow_difference = fd;
	down = true;
}

FloatPointPair Bristle::GetOffset()
{
	FloatPointPair ret = offset;
	ret.x += wander.x;
	ret.y += wander.y;
	return ret;
}

FloatPointPair Bristle::GetUnadjustedOffset()
{
	return offset;
}

float Bristle::GetFlowDiff()
{
	return flow_difference;
}

bool Bristle::AdjustOffset(FloatPointPair o)
{
	offset = o;
	return true;
}

bool Bristle::AdjustWander(FloatPointPair w)
{
	wander.x += w.x;
	wander.y += w.y;
	//	if ((wander.x * wander.x + wander.y + wander.y) > 2.0)
	//	{
	wander.x = (float)0.9 * wander.x;
	wander.y = (float)0.9 * wander.y;
	//	}
	return true;
}

bool Bristle::SetLast(FloatPointPair loc)
{
	last_loc = loc;
	return true;
}

FloatPointPair Bristle::GetLast()
{
	return last_loc;
}

bool Bristle::GetBristleDown()
{
	return down;
}

bool Bristle::SetBristleDown(bool d)
{
	down = d;
	return true;
}

Brush::Brush(FloatPointPair start, Color c, Color sec, float w, float d, Paint_Properties prop, Paper* paper, int pigment_index)
{
	paint_prop = prop;
	location = start;
	color = c;
	second = sec;
	bristle_kernel = NULL;
	watercolor = prop.watercolor;
	watercolor_paper = paper;
	watercolor_pigment_index = pigment_index;

	if (watercolor)
	{
		if (NULL == watercolor_paper)
		{
			throw std::runtime_error("Called Brush constructor for watercolor with NULL paper pointer.\n");
			return;
		}
		if (pigment_index < 0)
		{
			Pigment* pgmnt = new Pigment(watercolor_paper->GetWidth(), watercolor_paper->GetHeight(), &color, 0.2, 0.01, 1.0, 0.31);
			watercolor_pigment_index = watercolor_paper->SetPigment(pgmnt);
		}
	}

	if (paint_prop.brush_width_override)
	{
		brush_width = paint_prop.brush_width;
		num_bristles = brush_width * paint_prop.bristles;
	}
	else {
		brush_width = w * paint_prop.brush_width_factor;
		num_bristles = w * paint_prop.bristles;
	}
	brush_depth = d;

	shape = paint_prop.shape;
	orientation = 0.0;
	float flow_scale;

	if (paint_prop.paint_scale <= 0)
	{
		paint_prop.paint_scale = 0.01f;
	}
	if (paint_prop.flow <= 0)
	{
		paint_prop.flow = 1.0;
	}
	if (paint_prop.flow > 100)
	{
		paint_prop.flow = 100.0;
	}

	flow_scale = 1.0;


	paint_prop.flow = paint_prop.flow * flow_scale;
	paint_prop.flow_variation = paint_prop.flow_variation * flow_scale;
	if ((paint_prop.flow_variation + paint_prop.flow) > 100)
	{
		paint_prop.flow_variation = 100 - paint_prop.flow;
	}
	if (paint_prop.flow_variation < 0)
	{
		paint_prop.flow_variation = 0;
	}

	float radius = 0;
	float direction = 0;
	float fd = 0;
	std::random_device rd;
	std::mt19937 gen(rd());

	std::uniform_real_distribution<> rad_dist(0, 1.0);
	std::uniform_real_distribution<> direction_dist(0, 6.2832);
	std::uniform_real_distribution<> flow_dist(-paint_prop.flow_variation, paint_prop.flow_variation);
	std::uniform_real_distribution<> wide_dist(-brush_width, brush_width);
	std::uniform_real_distribution<> deep_dist(-brush_depth, brush_depth);

	bristle_kernel = (float**)malloc(sizeof(float*) * 3);
	if (NULL == bristle_kernel)
	{
		throw std::runtime_error("Error allocating memory for bristle_kernel.\n");
		return;
	}
	for (int i = 0; i <= 2; ++i)
	{
		bristle_kernel[i] = (float*)malloc(sizeof(float) * 3);
		if (NULL == bristle_kernel[i])
		{
			throw std::runtime_error("Error allocating memory for bristle_kernel.\n");
			return;
		}
		for (int j = 0; j <= 2; ++j)
		{
			bristle_kernel[i][j] = exp(-sqrt(i * i + j * j) * prop.bristle_thin_factor);
		}
	}

	bristles.clear();
	FloatPointPair o = { 0, 0 };
	if (shape_test == shape)
	{
		num_bristles = 2 * brush_width;
	}
	for (int i = 0; i < num_bristles; ++i)
	{
		if (shape_round == shape)
		{
			radius = brush_width * sqrt(rad_dist(gen));
			direction = direction_dist(gen);
			fd = flow_dist(gen);
			o.x = radius * cos(direction);
			o.y = radius * sin(direction);
		}
		else if (shape_straight == shape)
		{
			fd = flow_dist(gen);
			o.x = deep_dist(gen);
			o.y = wide_dist(gen);
		}
		else if (shape_test == shape)
		{
			fd = 0;
			o.x = 0;
			o.y = -brush_width + 2 * brush_width * ((float)i / (float)num_bristles);
		}
		else
		{
			num_bristles = 1;
			fd = 0;
			o.x = 0;
			o.y = 0;
		}

		bristles.push_back(new Bristle(o, fd));
		if (NULL == bristles.back())
		{
			throw std::runtime_error("Failed to allocate memory for bristle.\n");
			return;
		}
	}
}

Brush::~Brush()
{
	if (NULL != bristle_kernel)
	{
		for (int i = 0; i <= 2; ++i)
		{
			if (NULL != bristle_kernel[i])
			{
				free(bristle_kernel[i]);
			}
		}
		free(bristle_kernel);
	}
	if (bristles.size() > 0)
	{
		for (std::vector<Bristle*>::iterator it = bristles.begin(); it != bristles.end(); ++it)
		{
			if (NULL != *it)
			{
				delete (*it);
			}
		}
		bristles.clear();
	}
}


bool Brush::MoveTo(FloatPointPair loc)
{
	//location.x = loc.x * paint_scale;
	//location.y = loc.y * paint_scale;
	location = loc;
	return true;
}

bool Brush::PaintTo2(FloatPointPair p2, FloatPointPair o2, float* data, int width, int height, SPixelData* mask, int mask_value, bool use_mask, float rad1, float rad2, bool begin, SPixelData* extinguish_mask)
{
	FloatPointPair p1 = { 0, 0 }, p_temp = { 0, 0 }, direction = { 0, 0 };
	float gradient;
	float dr;
	float r1 = rad1;
	float r2 = rad2;
	float curve_adjustment = paint_prop.paint_scale;
	float dorient = 0.0;
	float orient1 = orientation;
	float orient2;
	float delta_orientation;

	if (false == paint_prop.radius_variation)
	{
		r1 = -1;
		r2 = -1;
	}
	orient2 = CalculateOrientation(o2);

	delta_orientation = OrientationDifference(orient1, orient2);
	if (paint_prop.glitch1)
	{
		curve_adjustment = paint_prop.paint_scale / 4.0f;
	}


	if (rad1 < 0)
	{
		if (rad2 < 0)
		{
			r1 = brush_width;
		}
		else {
			r1 = rad2;
		}
	}
	if (rad2 < 0)
	{
		if (rad1 < 0)
		{
			r2 = brush_width;
		}
		else {
			r2 = rad1;
		}
	}
	int iterations;
	FloatPointPair current_location;
	p1 = location;
	if ((p1.x == p2.x) && (p1.y == p2.y))
	{
		if (use_mask)
		{
			Dab3(o2, data, width, height, mask, mask_value, r1, begin, extinguish_mask);
		}
		else {
			Dab3(o2, data, width, height, NULL, mask_value, r1, begin, extinguish_mask);
		}
		return true;
	}
	float dx, dy;
	dx = p2.x - p1.x;
	dy = p2.y - p1.y;
	float dist = sqrt(dx * dx + dy * dy);
	float factor;
	if (r2 > r1)
	{
		factor = 1.0f + r2 * abs(tan(abs(delta_orientation))) / dist;
	}
	else {
		factor = 1.0f + r1 * abs(tan(abs(delta_orientation))) / dist;
	}
	if (factor < 1.0f)
	{
		factor = 1.0f;
	}
	else if (factor > MAX_CURVE_FACTOR)
	{
		factor = MAX_CURVE_FACTOR;
	}

	curve_adjustment = curve_adjustment * factor;
	float inverse_curve_adjustment = 1.0f / curve_adjustment;

	bool steep = (abs(p1.y - p2.y) > abs(p1.x - p2.x));
	int step;
	if (steep)
	{
		if (p1.y > p2.y)
		{
			step = -1;
		}
		else {
			step = 1;
		}
		gradient = step * (float)(p2.x - p1.x) / (float)((p2.y - p1.y) * curve_adjustment);
		dr = step * (float)(r2 - r1) / (float)((p2.y - p1.y) * curve_adjustment);
		dorient = step * delta_orientation / (float)((p2.y - p1.y) * curve_adjustment);
		iterations = abs(p2.y - p1.y) * curve_adjustment;
	}
	else {
		if (p1.x > p2.x)
		{
			step = -1;
		}
		else {
			step = 1;
		}
		gradient = step * (float)(p2.y - p1.y) / (float)((p2.x - p1.x) * curve_adjustment);
		dr = step * (float)(r2 - r1) / (float)((p2.x - p1.x) * curve_adjustment);
		dorient = step * delta_orientation / (float)((p2.x - p1.x) * curve_adjustment);
		iterations = abs(p2.x - p1.x) * curve_adjustment;
	}
	current_location = p1;

	direction.x = cos(orient1);
	direction.y = sin(orient1);
	if (!use_mask)
	{
		mask = NULL;
		//		mask_value = 0;
	}

	Dab3(direction, data, width, height, mask, mask_value, r1, begin, extinguish_mask);
	for (int i = 0; i < iterations; ++i)
	{
		if (steep)
		{
			current_location.y = current_location.y + step * inverse_curve_adjustment;
			current_location.x = current_location.x + gradient;
		}
		else {
			current_location.x = current_location.x + step * inverse_curve_adjustment;
			current_location.y = current_location.y + gradient;
		}
		r1 += dr;
		orient1 += dorient;
		if (orient1 > 2 * M_PI)
		{
			orient1 -= 2 * M_PI;
		}
		else if (orient1 < 0)
		{
			orient1 += 2 * M_PI;
		}
		direction.x = cos(orient1);
		direction.y = sin(orient1);
		MoveTo(current_location);
		SetOrientation(orient1);
		Dab3(direction, data, width, height, mask, mask_value, r1, false, extinguish_mask);
	}
	MoveTo(p2);
	SetOrientation(orient2);
	return true;
}

bool Brush::ChangeColor(Color c, Color sec)
{
	color = c;
	second = sec;
	return true;
}

bool Brush::Dab3(FloatPointPair direction, float* data, int width, int height, SPixelData* mask, int mask_value, float spot_radius, bool begin, SPixelData* extinguish_mask)
{
	// Use direction to measure whether the movement has been forward.
	if ((abs(direction.x) > EFFECTIVE_ZERO) && (abs(direction.y) > EFFECTIVE_ZERO))
	{
		float magnitude = sqrt(direction.x * direction.x + direction.y * direction.y);
		direction.x = direction.x / magnitude;
		direction.y = direction.y / magnitude;
	}

	float scaled_spot_radius_squared = spot_radius * spot_radius * paint_prop.paint_scale * paint_prop.paint_scale;
	bool new_stroke = false;
	bool extinguish = (NULL != extinguish_mask);

	std::random_device rd;
	std::mt19937 gen(rd());

	std::uniform_real_distribution<> Adjust(0, 1);
	std::uniform_real_distribution<> Adjust_x(-1.0, 1.0);
	std::uniform_real_distribution<> Adjust_y(-1.0, 1.0);

	Bristle* bristle;
	FloatPointPair oriented_location = { 0, 0 };  // Location of bristle based on bristle offset and orientation of brush (location relative to brush center).

	float c_value = cos(orientation);
	float s_value = sin(orientation);
	int mask_width = 0;
	int mask_height = 0;
	std::array<FloatPointPair, 4> points = { location, location, location, location }; // Points used for extinguishing on the extinguish_mask.
	std::array<float, 4> corner_dist = { 0.0, 0.0, 0.0, 0.0 };  // Distance towards each corner of points.
	if (NULL != mask)
	{
		mask_width = mask->GetWidth();
		mask_height = mask->GetHeight();
	}

	int count = 0;
	for (std::vector<Bristle*>::iterator it = bristles.begin(); it != bristles.end(); ++it)
	{
		count++;
		new_stroke = false; // Default.
		Color bristle_color;
		if (0 == count % 50)
		{
			bristle_color = second;
		}
		else {
			bristle_color = color;
		}
		bristle = *it;
		FloatPointPair offset = bristle->GetOffset();
		oriented_location.x = paint_prop.paint_scale * (c_value * offset.x - s_value * offset.y);
		oriented_location.y = paint_prop.paint_scale * (c_value * offset.y + s_value * offset.x);
		FloatPointPair raw_offset = bristle->GetUnadjustedOffset();  // Offset location of bristle without taking into account wander (relative to brush center).
		float x = c_value * raw_offset.x - s_value * raw_offset.y + location.x;  // Absolute location of bristle (not relative to brush center), without wander.
		float y = c_value * raw_offset.y + s_value * raw_offset.x + location.y;
		FloatPointPair bristle_location = { 0, 0 };  // Absolute location of bristle (not relative to brush center), with wander.
		bristle_location.x = paint_prop.paint_scale * location.x + oriented_location.x;
		bristle_location.y = paint_prop.paint_scale * location.y + oriented_location.y;
		FloatPointPair last;  // Where was the last place the bristle was painted?
		float movement_distance = 0.0f;  // How far has bristle moved since last painted spot?
		float dot_product = 0.0f; // What is the direction of bristle movement when compared to brush movement?  Positive is same general direction.
		float flow_diff = bristle->GetFlowDiff();
		float channels[3] = { 0.0f, 0.0f, 0.0f };
		if (begin)
		{
			new_stroke = true;
			bristle->SetLast(bristle_location);
			bristle->SetBristleDown(true);
		}
		else
		{
			last = bristle->GetLast();
			FloatPointPair bristle_movement = { 0, 0 };
			bristle_movement.x = bristle_location.x - last.x;
			bristle_movement.y = bristle_location.y - last.y;
			movement_distance = bristle_movement.x * bristle_movement.x + bristle_movement.y * bristle_movement.y;
			dot_product = bristle_movement.x * direction.x + bristle_movement.y * direction.y;
		}

		if (((NULL == mask) || ((x >= 0) && (x < mask_width) && (y >= 0) && (y < mask_height) && (mask_value == mask->GetPixel(x, y)))) &&
			((spot_radius < 0) || ((oriented_location.x * oriented_location.x + oriented_location.y * oriented_location.y) < scaled_spot_radius_squared)) &&
			(dot_product >= 0.0f))
		{
			if (false == bristle->GetBristleDown())  // Newly in writeable area.
			{
				new_stroke = true;
				bristle->SetLast(bristle_location);
				bristle->SetBristleDown(true);
			}
		}
		else {  // Outside of mask or writeable area.
			bristle->SetBristleDown(false);
		}
		if ((bristle->GetBristleDown()) &&  // First, make sure bristle is down.
			(new_stroke || ((movement_distance >= 1.0))))  // But we also need to either have a new stroke or sufficient movement distance.
		{
			bristle->SetLast(bristle_location);
			if (watercolor)
			{
//				watercolor_paper->Dab((int)x, (int)y, 5, 0.6, (flow_diff + paint_prop.flow) / 5000.0, watercolor_pigment_index);  // *** Try doing a single pass of the Superpixel with uniform concentration, and velocity from path. ***
				// Just need to adjust velocity.
				watercolor_paper->SetVelocity((int)x, (int)y, direction.x/1000.0f, direction.y/1000.0f, true);
			}
			else {
				for (int i = -2; i <= 2; ++i)
				{
					for (int j = -2; j <= 2; ++j)
					{
						float adjustment;
						x = bristle_location.x + i;
						y = bristle_location.y + j;
						if (paint_prop.sub_pixel)
						{
							adjustment = KernelAdjustment(i, j, x, y);
						}
						else {
							adjustment = bristle_kernel[abs(i)][abs(j)];
						}

						float adjusted_flow = (flow_diff + paint_prop.flow) * adjustment;
						if ((x >= 0) && (x < width) && (y >= 0) && (y < height))
						{

							for (int color_index = 0; color_index < 3; ++color_index)
							{
								long pos = (int)y * width * 3 + (int)x * 3 + color_index;
								if (data[pos] > bristle_color.channel[color_index]) // Brighter background
								{
									channels[color_index] = (adjusted_flow * (float)bristle_color.channel[color_index] + (100.0f - adjusted_flow) * (float)data[pos]) / 100.0f;
									if (channels[color_index] < bristle_color.channel[color_index])
									{
										data[pos] = bristle_color.channel[color_index];
									}
									else if (channels[color_index] > 255)
									{
										data[pos] = 255;
									}
									else {
										if (paint_prop.glitch2)
										{
											data[pos] = (float)((unsigned char)channels[color_index]);
										}
										else {
											data[pos] = channels[color_index];
										}
									}
								}
								else { // Dimmer background
									channels[color_index] = (adjusted_flow * (float)bristle_color.channel[color_index] + (100.0f - adjusted_flow) * (float)data[pos]) / 100.0f;
									if (channels[color_index] < 0)
									{
										data[pos] = 0;
									}
									else if (channels[color_index] > bristle_color.channel[color_index])
									{
										data[pos] = bristle_color.channel[color_index];
									}
									else {
										if (paint_prop.glitch2)
										{
											data[pos] = (float)((unsigned char)channels[color_index]);
										}
										else {
											data[pos] = channels[color_index];
										}
									}
								}
							}

						}
					}
				}
			}
		}
		if (Adjust(gen) < 0.1)
		{
			offset.x = Adjust_x(gen);
			offset.y = Adjust_y(gen);
			bristle->AdjustWander(offset);
		}
	}
	if (extinguish)
	{
		for (int i = 0; i < 4; i++)
		{
			PointPair corner_location;
			float mod_brush_width = brush_width - 2.5;
			if (mod_brush_width < 5.0)
			{
				if (brush_width > 5.0)
				{
					mod_brush_width = 5.0;
				}
				else {
					mod_brush_width = brush_width;
				}
			}
			corner_location.x = -1 + (i & 2);
			corner_location.y = -1 + ((i + 1) & 2);
			points[i].x = corner_location.x * c_value * brush_depth - corner_location.y * s_value * mod_brush_width + location.x;
			points[i].y = corner_location.y * c_value * mod_brush_width + corner_location.x * s_value * brush_depth + location.y;
		}
		return ExtinguishQuadrilateral(points, width, height, extinguish_mask, mask_value);
	}
	return true;
}

bool Brush::PaintCorner2(Corner corner, float* data, int width, int height, SPixelData* mask, int mask_value, bool use_mask, SPixelData* extinguish_mask)
{
	float chord_length = sqrt((corner.p0.x - corner.p1.x) * (corner.p0.x - corner.p1.x) + (corner.p0.y - corner.p1.y) * (corner.p0.y - corner.p1.y));
	float leg1 = sqrt((corner.p0.x - corner.c0.x) * (corner.p0.x - corner.c0.x) + (corner.p0.y - corner.c0.y) * (corner.p0.y - corner.c0.y));
	float leg2 = sqrt((corner.c1.x - corner.c0.x) * (corner.c1.x - corner.c0.x) + (corner.c1.y - corner.c0.y) * (corner.c1.y - corner.c0.y));
	float leg3 = sqrt((corner.p1.x - corner.c1.x) * (corner.p1.x - corner.c1.x) + (corner.p1.y - corner.c1.y) * (corner.p1.y - corner.c1.y));
	float length_est = paint_prop.paint_scale * (chord_length + leg1 + leg2 + leg3) / 2.0f;
	float local_radius = -1;
	float dr = 0;
	FloatPointPair L1 = { 0, 0 }, L2 = { 0, 0 }, L3 = { 0, 0 }, M1 = { 0, 0 }, M2 = { 0, 0 }, N1 = { 0, 0 };
	FloatPointPair direction = { 0, 0 };

	int n = (int)length_est + 1;
	if (paint_prop.radius_variation)
	{
		local_radius = corner.radius_p0;
		dr = (corner.radius_p1 - corner.radius_p0) / (float)n;
	}
	MoveTo(corner.p0);
	direction.x = corner.c0.x - corner.p0.x;
	direction.y = corner.c0.y - corner.p0.y;
	if ((abs(direction.x) <= EFFECTIVE_ZERO) && (abs(direction.y) <= EFFECTIVE_ZERO))
	{
		direction.x = corner.c1.x - corner.p0.x;
		direction.y = corner.c1.y - corner.p0.y;
	}
	SetOrientation(direction);

	for (int i = 0; i < n; ++i)
	{
		float p = (float)i / (float)n; // Parameter goes from 0 to 1.0;
		float np = 1.0f - p;
		L1.x = corner.c0.x * p + corner.p0.x * np;
		L1.y = corner.c0.y * p + corner.p0.y * np;
		L2.x = corner.c1.x * p + corner.c0.x * np;
		L2.y = corner.c1.y * p + corner.c0.y * np;
		L3.x = corner.p1.x * p + corner.c1.x * np;
		L3.y = corner.p1.y * p + corner.c1.y * np;
		M1.x = L2.x * p + L1.x * np;
		M1.y = L2.y * p + L1.y * np;
		M2.x = L3.x * p + L2.x * np;
		M2.y = L3.y * p + L2.y * np;
		N1.x = M2.x * p + M1.x * np;
		N1.y = M2.y * p + M1.y * np;
		direction.x = M2.x - M1.x;
		direction.y = M2.y - M1.y;
		PaintTo2(N1, direction, data, width, height, mask, mask_value, use_mask, local_radius, (local_radius + dr), 0 == i, extinguish_mask);

		local_radius += dr;
	}
	return true;
}

bool Brush::ExtinguishCorner(Corner corner, float* data, int width, int height, SPixelData* mask, int mask_value)
{
	// *** Probably don't need now.

	float chord_length = sqrt((corner.p0.x - corner.p1.x) * (corner.p0.x - corner.p1.x) + (corner.p0.y - corner.p1.y) * (corner.p0.y - corner.p1.y));
	float leg1 = sqrt((corner.p0.x - corner.c0.x) * (corner.p0.x - corner.c0.x) + (corner.p0.y - corner.c0.y) * (corner.p0.y - corner.c0.y));
	float leg2 = sqrt((corner.c1.x - corner.c0.x) * (corner.c1.x - corner.c0.x) + (corner.c1.y - corner.c0.y) * (corner.c1.y - corner.c0.y));
	float leg3 = sqrt((corner.p1.x - corner.c1.x) * (corner.p1.x - corner.c1.x) + (corner.p1.y - corner.c1.y) * (corner.p1.y - corner.c1.y));
	float length_est = paint_prop.paint_scale * (chord_length + leg1 + leg2 + leg3) / 2.0f;
	float local_radius = -1;
	float dr = 0;
	FloatPointPair L1 = { 0, 0 }, L2 = { 0, 0 }, L3 = { 0, 0 }, M1 = { 0, 0 }, M2 = { 0, 0 }, N1 = { 0, 0 };
	FloatPointPair direction = { 0, 0 };
	//float flow_adjustment = 1.0 / paint_scale;

	int n = (int)length_est + 1;
	if (paint_prop.radius_variation)
	{
		local_radius = corner.radius_p0;
		dr = (corner.radius_p1 - corner.radius_p0) / (float)n;
	}
	MoveTo(corner.p0);
	direction.x = corner.c0.x - corner.p0.x;
	direction.y = corner.c0.y - corner.p0.y;
	if ((abs(direction.x) <= EFFECTIVE_ZERO) && (abs(direction.y) <= EFFECTIVE_ZERO))
	{
		direction.x = corner.c1.x - corner.p0.x;
		direction.y = corner.c1.y - corner.p0.y;
	}
	SetOrientation(direction);

	for (int i = 0; i < n; ++i)
	{
		float p = (float)i / (float)n; // Parameter goes from 0 to 1.0;
		float np = 1.0f - p;
		L1.x = corner.c0.x * p + corner.p0.x * np;
		L1.y = corner.c0.y * p + corner.p0.y * np;
		L2.x = corner.c1.x * p + corner.c0.x * np;
		L2.y = corner.c1.y * p + corner.c0.y * np;
		L3.x = corner.p1.x * p + corner.c1.x * np;
		L3.y = corner.p1.y * p + corner.c1.y * np;
		M1.x = L2.x * p + L1.x * np;
		M1.y = L2.y * p + L1.y * np;
		M2.x = L3.x * p + L2.x * np;
		M2.y = L3.y * p + L2.y * np;
		N1.x = M2.x * p + M1.x * np;
		N1.y = M2.y * p + M1.y * np;
		direction.x = M2.x - M1.x;
		direction.y = M2.y - M1.y;
		// *** Figure out how to sweep out this.  Maybe create new extinguish version of PaintTo2. ***
//		PaintTo2(N1, direction, data, width, height, mask, mask_value, use_mask, local_radius, (local_radius + dr), 0 == i);



		local_radius += dr;
	}
	return true;
}

bool Brush::SetOrientation(FloatPointPair o)
{

	if (abs(o.x) > EFFECTIVE_ZERO)
	{
		SetOrientation(atan2(o.y, o.x));
	}
	else {
		if (o.y > 0)
		{
			orientation = M_PI / 2;
		}
		else {
			orientation = 3 * M_PI / 2;
		}
	}
	return true;
}

bool Brush::SetOrientation(float o)
{
	while (o > 2 * M_PI)
	{
		o -= 2 * M_PI;
	}
	while (o < 0)
	{
		o += 2 * M_PI;
	}
	orientation = o;
	return true;
}

float Brush::GetOrientation()
{
	return orientation;
}

float Brush::CalculateOrientation(FloatPointPair direction)
{
	// Taking a cartesian direction vector as input, calculate an orientation in radians (where zero points to 1.0,0.0 and positive values from there move in the positive y direction).
	// The range of outputs are 0 to 2*M_PI.
	float ret = 0.0f;

	if (abs(direction.x) > EFFECTIVE_ZERO)
	{
		ret = atan2(direction.y, direction.x);
		if (ret < 0)
		{
			ret += 2 * M_PI;
		}
	}
	else {
		if (direction.y > 0)
		{
			ret = 1.5708f;
		}
		else {
			ret = 4.7124f;
		}
	}
	return ret;
}

float Brush::OrientationDifference(float o1, float o2)
{
	// Unlike the orientation variable, this return value goes from -M_PI to +M_PI.
	float ret = 0.0f;
	ret = o2 - o1;
	if (abs(ret) > M_PI)
	{
		if (ret < 0)
		{
			ret += 2 * M_PI;
		}
		else {
			ret -= 2 * M_PI;
		}
	}
	return ret;
}

float Brush::KernelAdjustment(int i, int j, float x, float y)
{
	float ret = 0;
	float mx = x - (int)x; // Mod 1 of the x centerpoint.
	float my = y - (int)y; // Mod 1 of the y centerpoint.
	float nmx = 1.0f - mx;
	float nmy = 1.0f - my;

	float f00, f10, f01, f11;
	f00 = nmx * nmy;
	f01 = nmx * my;
	f10 = mx * nmy;
	f11 = mx * my;

	int k = i;
	int l = j;
	if ((k > -3) && (k < 3) && (l > -3) && (l < 3))
	{
		ret += bristle_kernel[abs(k)][abs(l)] * f00;
	}
	k--;
	if ((k > -3) && (k < 3) && (l > -3) && (l < 3))
	{
		ret += bristle_kernel[abs(k)][abs(l)] * f10;
	}
	l--;
	if ((k > -3) && (k < 3) && (l > -3) && (l < 3))
	{
		ret += bristle_kernel[abs(k)][abs(l)] * f11;
	}
	k++;
	if ((k > -3) && (k < 3) && (l > -3) && (l < 3))
	{
		ret += bristle_kernel[abs(k)][abs(l)] * f01;
	}
	return ret;
}

bool Brush::ExtinguishQuadrilateral(std::array<FloatPointPair, 4> points, int width, int height, SPixelData* extinguish_mask, int value)
{
	// Takes a convex quadrilateral, defined by sequential "points".  For each point on extinguish_mask that is in the interior of the polygon and with a value of "value", set value to zero.

	int x0, y0, x1, y1;  // Bounds of the quadrilateral.
	std::array<float, 4> m = { 0.0, 0.0, 0.0, 0.0 }, b = { 0.0, 0.0, 0.0, 0.0 }; // The slope and intercept values for the four line segments.  Here, the line equation is: y = m*x + b.
	std::array<bool, 4> steep = { false, false, false, false }; // True for line segments that are more vertical than horizontal.
	std::array<bool, 4> direction = { true, true, true, true }; // True indicates that interior is in the positive y direction (if not steep) or the positive x direction (if steep).
	// First, extinguish along all line segments.
	for (int i = 0; i < 4; i++)
	{
		int j = (i + 1) % 4;
		std::array<FloatPointPair, 2> ls_points;
		ls_points[0] = points[i];
		ls_points[1] = points[j];
		if (!ExtinguishLineSegment(ls_points, width, height, extinguish_mask, value))
		{
			throw std::runtime_error("Error calling ExtinguishLineSegment.\n");
			return false;
		}
	}
	// Next, find the lowest values of x and y.
	x0 = width;
	y0 = height;
	for (int i = 0; i < 4; i++)
	{
		if (points[i].x < x0)
		{
			x0 = points[i].x;
		}
		if (points[i].y < y0)
		{
			y0 = points[i].y;
		}
	}
	// Next, find largest values of x and y.
	x1 = 0;
	y1 = 0;
	for (int i = 0; i < 4; i++)
	{
		if (points[i].x > x1)
		{
			x1 = points[i].x;
		}
		if (points[i].y > y1)
		{
			y1 = points[i].y;
		}
	}
	// Calculate line properties for each line segment.
	for (int i = 0; i < 4; i++)
	{
		int j = (i + 1) % 4; // j is the index for the subsequent vertex, the other side of the line segment that starts with i.
		if (abs(points[j].x - points[i].x) < abs(points[j].y - points[i].y))
		{
			float denominator = points[j].y - points[i].y;
			if (abs(denominator) < EFFECTIVE_ZERO)
			{
				// Points i and j are the same.  Move point j to halfway between it and the next point.  If the next point is the same, then there are no interior points.
				int k = (i + 2) % 4;
				if ((abs(points[k].x - points[j].x) < EFFECTIVE_ZERO) && (abs(points[k].y - points[j].y) < EFFECTIVE_ZERO))
				{
					return true;  // No error, just no interior points.
				}
				else {
					points[j].x = (points[j].x + points[k].x) / 2.0;
					points[j].y = (points[j].y + points[k].y) / 2.0;
					return ExtinguishQuadrilateral(points, width, height, extinguish_mask, value);
				}
			}
			steep[i] = true;
			m[i] = (points[j].x - points[i].x) / denominator;
			b[i] = points[i].x - m[i] * points[i].y;
		}
		else
		{
			float denominator = points[j].x - points[i].x;
			if (abs(denominator) < EFFECTIVE_ZERO)
			{
				// Points i and j are the same.  Move point j to halfway between it and the next point.  If the next point is the same, then there are no interior points.
				int k = (i + 2) % 4;
				if ((abs(points[k].x - points[j].x) < EFFECTIVE_ZERO) && (abs(points[k].y - points[j].y) < EFFECTIVE_ZERO))
				{
					return true;  // No error, just no interior points.
				}
				else {
					points[j].x = (points[j].x + points[k].x) / 2.0;
					points[j].y = (points[j].y + points[k].y) / 2.0;
					return ExtinguishQuadrilateral(points, width, height, extinguish_mask, value);
				}
			}
			steep[i] = false;
			m[i] = (points[j].y - points[i].y) / denominator;
			b[i] = points[i].y - m[i] * points[i].x;
		}
	}
	// Calculate directions (and validate concave structure of polygon).
	for (int i = 0; i < 4; i++)
	{
		int j = (i + 2) % 4; // This is a point that is not part of the line segment.
		int k = (i + 3) % 4; // This is the other point that is not part of the line segment.
		float j_dir, k_dir;
		if (steep[i])
		{
			j_dir = points[j].x - m[i] * points[j].y - b[i];
			k_dir = points[k].x - m[i] * points[k].y - b[i];
		}
		else {
			j_dir = points[j].y - m[i] * points[j].x - b[i];
			k_dir = points[k].y - m[i] * points[k].x - b[i];
		}
		if ((j_dir >= 0.0) && (k_dir >= 0.0))
		{
			direction[i] = true;
		}
		else if ((j_dir <= 0.0) && (k_dir <= 0.0))
		{
			direction[i] = false;
		}
		else {  // *** Here. ***
			throw std::runtime_error("Non-convex polygon supplied to ExtinguishQuadrilateral.\n");
			return false;
		}
	}
	// Now step through all y values.
	for (int y = y0; y <= y1; y++)
	{
		for (int x = x0; x <= x1; x++)
		{
			if (extinguish_mask->GetPixel(x, y) == value)
			{
				bool interior = true;
				for (int i = 0; interior && (i < 4); i++)
				{
					if (steep[i])
					{
						interior = (((x - (m[i] * y + b[i])) >= 0) == direction[i]); // Compare the calculated direction of the point x,y with the interior direction from line segment i.
					}
					else {
						interior = (((y - (m[i] * x + b[i])) >= 0) == direction[i]); // Compare the calculated direction of the point x,y with the interior direction from line segment i.
					}
				}
				if (interior) // Point is in the interior.
				{
					extinguish_mask->SetPixel(x, y, 0); // Extinguish this pixel from the set.
				}
			}
		}
	}
	return true;
}

bool Brush::ExtinguishLineSegment(std::array<FloatPointPair, 2> points, int width, int height, SPixelData* extinguish_mask, int value)
{
	// Two "points" define a line segment.  Set all pixels with value of "value" in the line between them on the "extinguish_mask" to zero.
	float x, y, dx, dy, dist;
	// First, extinguish the vertices.
	for (int i = 0; i < 2; i++)
	{
		if (extinguish_mask->GetPixel(points[i].x, points[i].y) == value)
		{
			extinguish_mask->SetPixel(points[i].x, points[i].y, 0);
		}
	}

	dx = points[1].x - points[0].x;
	dy = points[1].y - points[0].y;
	dist = sqrt(dx * dx + dy * dy);
	int iterations = (int)dist + 1;
	dx = dx / (float)iterations;
	dy = dy / (float)iterations;
	x = points[0].x;
	y = points[0].y;
	for (int i = 0; i < iterations; i++)
	{
		x += dx;
		y += dy;
		if (extinguish_mask->GetPixel(x, y) == value)
		{
			extinguish_mask->SetPixel(x, y, 0);
		}
	}
	return true;
}
