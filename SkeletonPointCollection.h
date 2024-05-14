#pragma once

// Copyright (c) 2023-2024 Steve Young
// Licensed under the MIT License

#include <vector>
#include <set>

class SuperPixel;
class SPixelData;
struct PointPair;

struct Skeleton_Link
{
	int end1; // One endpoint or intersection point.
	int end2; // The other endpoint or intersection point.
	std::set<int> link_points; // All connecting link points adjacent to end1 or end2. (Not all points, just those around end1 or end2.)
	std::vector<int> line; // sequential points from end1 to end2
	float distance;
	bool used_in_long_path = false;
	bool operator< (const Skeleton_Link& rhs) const
	{
		return end1 < rhs.end1;
	}
};

struct Skeleton_Point
{
	int point; // The point that shows up in endpoint_set or intersection_set.
	bool endpoint;  // True if an endpoint, false if an intersection.
	Skeleton_Link* dir[8]; // dir[0] is above, each incremental index is around the pixel in the clockwise direction.
	float temp_dist; // Temporary distance value used for finding the shortest path between two Skeleton_Points.
	bool operator< (const Skeleton_Point& rhs) const
	{
		return point < rhs.point;
	}
};


struct SkeletonPointIndex
{
	int point;
	Skeleton_Point* sp;
	bool operator< (const SkeletonPointIndex& rhs) const
	{
		return point < rhs.point;
	}
};

class SkeletonPointCollection
{
private:
	SuperPixel* skeleton_superpixel;  // This is the SuperPixel that is thinned to a skeleton.  Not the full SuperPixel representation.
	std::set<SkeletonPointIndex> skeleton_point_set;
	int width;
	int height;
	SPixelData* pixdata;
	std::vector<int> longest_path;
	std::set<std::vector<int>> longest_paths;

public:
	SkeletonPointCollection(SuperPixel* sp);
	~SkeletonPointCollection();

	bool FindSkeletonPoints();
	bool FindSkeletonLinks();
	bool FindLongestPaths(bool just_one = false);
	bool RecalcDistances();
	int FindDirSmallestTempDist(Skeleton_Point* sp);
	bool AddSkeletonEndpoint(PointPair p);
	bool AddSkeletonIntersection(PointPair p);
	std::set<int> GetEndpoints();
	std::set<int> GetIntersections();
	std::set<int> GetAdjacentIntersections(int pos);
	bool TestEndpointOrIntersection(int pos);

	Skeleton_Point* GetSPByIdentifier(int key);
	bool InsertPointIntoIndex(Skeleton_Point* sp);
	int GetDirection(int p1, int p2, int width);
	std::vector<int>* GetLongestPath();
	std::set<std::vector<int>> GetLongestPaths();

	bool ClearPoints();
	bool ResetTempDists();
	bool ClearLongPathFlags();
};
