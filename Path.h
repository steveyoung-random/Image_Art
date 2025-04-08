#pragma once

// Copyright (c) 2023-2025 Steve Young
// Licensed under the MIT License

#include "SPixelData.h"
#include "GradData.h"
#include "SuperPixel.h"
#include "General.h"
#include "stb_image_write.h"
#include <vector>
#include <set>
#include <sstream>

#define LINE_DISTANCE 1.0
#define VERTEX_ADJUSTMENT
#define CORNER_ALPHA_MAX 1.3  // This seems to produce smoother outputs at 1.3 rather than the 1.0 of the original default.
#define CORNER_ALPHA_MIN 0.55
#define SMOOTH_BOX_SIZE 0.25  // This seems to produce smoother outputs at 0.25 than the 0.5 from the original default.

struct Segment {
	int j;
	float penalty;
};

struct Segments {
	std::vector<Segment> s;
	int vertex;
	int prev_vertex;
	int num_segments;
	float cumulative_penalty;
	bool meeting_point = false;
};

struct PolygonVertex {
	FloatPointPair point;
	bool meeting_point;
};

class Path {
private:
	std::vector<PointPair> points;
	std::set<int> point_set;
	std::vector<PointPair> fine_points;  // Rather than record sub-pixel location, record pixel location of pixel down and to the right of the point.
	std::vector<PointPair> path_points; // Either points or fine_points, depending on which is being used.  Set in CalcForwardSegments().
	Path* next = NULL;
	Path* prev = NULL;
	Path* subpath = NULL;
	RectQuad boundingbox;
	bool forward;
	std::string color_string;
	SuperPixel* superpixel = NULL;
	std::vector<Segments> potential_segs;
	std::vector<PolygonVertex> optimal_segments;
	std::vector<Corner> optimal_curve;
	bool closed;
	bool begin_endpoint;
	bool end_endpoint;

public:
	Path(SuperPixel* sp, PointPair start, bool use_meeting_points = true, bool create_polygons = true, bool fine = true, std::vector<int>* line = NULL, bool beg_ep = false, bool end_ep = false, bool glitch3 = false);
	~Path();
	std::vector<PointPair> GetPoints();
	std::set<int> GetPointSet();
	std::vector<PointPair> GetFinePoints();
	std::vector<PolygonVertex> GetOptSegments();
	std::vector<Corner> GetCurve();
	FloatPointPair AdjustVertex(int prev_point, int current_point, int next_point);
	RectQuad GetBox();
	bool GetForward();
	bool GetClosed();
	std::string GetColor();
	Path* GetByEdgePixel(PointPair point);
	Path* GetHead();
	Path* GetNext();
	Path* GetPrev();
	Path* GetTail();
	Path* GetSubpathHead();
	Path* GetSubpathTail();
	FloatPointPair LinearRegression(int point1, int point2);
	FloatPointPair SquareNeighborhood(FloatPointPair vertex, FloatPointPair point);

	bool CalcPoints(PointPair start, bool use_meeting_points = true);
	bool CalcForwardSegments(bool use_meeting_points = true, bool fine = true);
	bool CalcOptimalSegments();
	bool CalcSegmentPenalties();
	bool Curve(bool glitch3 = false);

	bool InsertNext(Path* n);
	bool InsertPrev(Path* p);
	bool InsertSubpath(Path* p);

	bool SetNext(Path* n);
	bool SetPrev(Path* p);

	bool MoveSuperPixel(SuperPixel* newSP);

	Corner CalculateCorner(FloatPointPair p0, FloatPointPair p1, FloatPointPair p2, float alpha_max = CORNER_ALPHA_MAX, bool extend_p0 = false, bool extend_p2 = false, bool glitch3 = false);
};

