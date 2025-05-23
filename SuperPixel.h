#pragma once

// Copyright (c) 2023-2025 Steve Young
// Licensed under the MIT License

#include "SPixelData.h"
#include "GradData.h"
#include "stb_image_write.h"
#include "Path.h"
#include <vector>
#include <set>
#include <sstream>
#include <algorithm>

class Path;
class WorkSpace;
class SkeletonPointCollection;

enum SuperPixelType{SPType_Plain, SPType_Processed, SPType_Skeleton};

struct Contrast_Histogram
{
	int identifier;
	unsigned int location;
	SuperPixel* sp;
	bool operator< (const Contrast_Histogram& rhs) const
	{
		return identifier < rhs.identifier;
	}
};

class SuperPixel {
private:
	SPixelData* pixeldata = NULL;
	GradData* gradientdata = NULL;
	SuperPixel* next;
	SuperPixel* prev;
	int identifier;
	int size, prevsize, level_complete;
	PointPair seed;
	FloatPointPair centroid;
	RectQuad boundingbox;
	std::set<int> EdgePixels;
	std::set<int> Neighbors;
	std::set<int> Vertices;
	Color AveColor;
	float AveError;
	int num_paths;
	Path* path_list_head = NULL;
	Path* path_list_tail = NULL;
	WorkSpace* workspace = NULL;
	std::set<int> skeleton_endpoint_set;
	std::set<int> skeleton_intersection_set;
	SkeletonPointCollection* skeleton_collection = NULL;
	int color_bucket;
	SuperPixelType type;
	bool EdgePixelsCurrent;
	int image_width;
	int image_height;
	std::string fill_image;

public:
	SuperPixel(int id, GradData* graddat, SPixelData* pixdat, PointPair point, SuperPixel* n, SuperPixel* p, WorkSpace* ws = NULL, SuperPixelType t = SPType_Plain);
	SuperPixel(const SuperPixel& tsp, GradData* graddat, SPixelData* pixdat, SuperPixel* n, SuperPixel* p, WorkSpace* ws = NULL, SuperPixelType t = SPType_Plain);
	SuperPixel(int id, SPixelData* pixdat, PointPair point, SuperPixel* n, SuperPixel* p, WorkSpace* ws = NULL, SuperPixelType t = SPType_Plain);
	~SuperPixel();

	bool Absorb(SuperPixel* i, bool remove = true, bool manage_neighbors = false);
	bool AddNeighbor(int ne);
	bool AddVertex(int v);
	bool CalculateSize();
	bool CheckVertex(int v, bool check_neighbors = false);
	float ColorDifference(Color c1, Color c2);
	bool FindEdgePixels();
	bool FindPaths(bool use_meeting_points = true, bool polygon = true, bool fine = true);
	bool GenerateContrastImage(SPixelData* pdata, int radius);
	Color CalculateContrastColor(Color opposing);
	SuperPixel* DuplicateSuperPixelSet(SPixelData* sd = NULL);
	bool UpdateWorkspaceSuperPixelSet(WorkSpace* ws);
	bool DeletePathList();

	Color GetAveColor();
	float GetAveError();
	SuperPixel* GetByIdentifier(int key);
	FloatPointPair GetCentroid();
	int GetColorBucket();
	std::set<int>* GetEdgePixels();
	std::string GetFillImage();
	GradData* GetGradient();
	SuperPixel* GetHead();
	int GetIdentifier();
	std::set<int>* GetNeighbors();
	SuperPixel* GetNext();
	Path* GetPathHead();
	Path* GetPathTail();
	SPixelData* GetPixelData();
	SuperPixel* GetPrev();
	int GetPrevSize();
	PointPair GetSeed();
	int GetSize();
	int GetSetSize();
	SuperPixel* GetTail();
	SuperPixelType GetType();
	RectQuad GetWindow();
	WorkSpace* GetWorkspace();


	int Grow(unsigned char value, bool limit, bool mode, RectQuad box, SPixelData* mask, int mask_value, bool diagonals = true);
	bool Reset();

	Color SetAveColor(ImageData* image);
	Color SetAveColor(Color c);
	bool SetColorBucket(int b);
	bool SetFillImage(std::string fill);
	bool SetSize(int s);
	bool SetIdentifierandBox(int id, RectQuad box);
	bool SetNeighbors();
	bool SetPathHead(Path* path);
	bool SetPathTail();
	bool SetNext(SuperPixel* n);
	bool SetPrev(SuperPixel* p);
	bool SetPrevSize();
	bool SetWindow(RectQuad w);
	bool SetWorkspace(WorkSpace* ws);

	PointPair Split(ImageData* image, bool diagonals = false);
	bool SeparateDiscontinuousSuperPixel();
	SuperPixel* Thin(GradData* graddat, SPixelData* pixdat, SuperPixel* prev, bool only_one, bool glitch3 = false);
	bool Thin_Subiteration(int n, bool global_changed);

	PointPair Pos2XY(int pos);
	int XY2Pos(PointPair xy);

	std::vector<int> write_data();
	int Count();
	bool RemoveNeighbor(int ne);

	bool FindSkeletonPoints();
	bool AddSkeletonEndpoint(PointPair p);
	bool AddSkeletonIntersection(PointPair p);
	std::set<int> GetSkeletonEndpoints();
	std::set<int> GetSkeletonIntersections();

	bool ClearSkeletonPoints();

	std::string SaveState();
};
