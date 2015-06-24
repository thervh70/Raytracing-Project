#include <numeric>
#include "KDtree.h"
#include "mesh.h"
#include "raytracing.h"

AccelTreeNode treeRoot;
int treeDepth = 0;
int treeNodes = 1;

// Set the tree accuracy (choose values like 10, 100, 1000).
// only lower this if building the tree is taking too much time.
float TREE_ACCURACY = 1000.0f;

// Build a KD tree and store the planes with the triangles in a vector of vectors
void buildKDtree()
{
	std::cout << std::endl;
	std::cout << "Building tree ..." << std::endl;

	float xMin, xMax, yMin, yMax, zMin, zMax;

	// Calculate the bounding box of the scene
	std::vector<Triangle>::const_iterator it = MyMesh.triangles.begin();
	xMin = MyMesh.vertices[(*it).v[0]].p[0];
	xMax = MyMesh.vertices[(*it).v[0]].p[0];
	yMin = MyMesh.vertices[(*it).v[0]].p[1];
	yMax = MyMesh.vertices[(*it).v[0]].p[1];
	zMin = MyMesh.vertices[(*it).v[0]].p[2];
	zMax = MyMesh.vertices[(*it).v[0]].p[2];

	// Get the smallest and largest x, y and z values of the triangles of the mesh
	while (it != MyMesh.triangles.end()) {
		for (int i = 0; i < 3; ++i)
		{
			// check x
			if (MyMesh.vertices[(*it).v[i]].p[0] < xMin)
				xMin = MyMesh.vertices[(*it).v[i]].p[0];
			else if (MyMesh.vertices[(*it).v[i]].p[0]> xMax)
				xMax = MyMesh.vertices[(*it).v[i]].p[0];

			// check y
			if (MyMesh.vertices[(*it).v[i]].p[1] < yMin)
				yMin = MyMesh.vertices[(*it).v[i]].p[1];
			else if (MyMesh.vertices[(*it).v[i]].p[1] > yMax)
				yMax = MyMesh.vertices[(*it).v[i]].p[1];

			// check z
			if (MyMesh.vertices[(*it).v[i]].p[2] < zMin)
				zMin = MyMesh.vertices[(*it).v[i]].p[2];
			else if (MyMesh.vertices[(*it).v[i]].p[2] > zMax)
				zMax = MyMesh.vertices[(*it).v[i]].p[2];
		}
		++it;
	}

	// Assign the found values to the kd tree (with a small margin)
	// !!! IMPORTANT: NEED BETTER SOLUTION THAN INCREASING MARGIN TO ALLOW CAMERA TO BE OUTSIDE OF THE ROOT !!!
	std::vector<int> indexlist(MyMesh.triangles.size());
	std::iota(std::begin(indexlist), std::end(indexlist), 0); // Fill with 0, 1, ..., 99.
	treeRoot = AccelTreeNode(indexlist, xMin - 5.1f, xMax + 5.1f, yMin - 5.1f, yMax + 5.1f, zMin - 5.1f, zMax + 5.1f);

	// Add the root node to it's own parent list.
	treeRoot.parentList.push_back(&treeRoot);

	// Split the main node into smaller nodes
	splitSpaces(treeRoot, 0);

	std::cout << "... done building tree" << std::endl;
	std::cout << "Tree depth: " << treeDepth << std::endl;
	std::cout << "Amount of nodes: " << treeNodes << std::endl;
	std::cout << "Amount of triangles: " << MyMesh.triangles.size() << std::endl;
	std::cout << std::endl;
}

void splitSpaces(AccelTreeNode& tree, const int axis) {

	if (axis > treeDepth)
		++treeDepth;


	// split the KDtreeCube into subspaces and recursevily recall on those subspaces

	// stop if there are 10 or less triangles in the current level
	if (tree.triangles.size() < 11)
	{
		tree.leftChild = nullptr;
		tree.rightChild = nullptr;
		return;
	}


	// stop if subspace axis which is being divided is smaller than 0.0001
	if ((axis % 3 == 0 && tree.xEnd - tree.xStart < 0.0001) ||
		(axis % 3 == 1 && tree.yEnd - tree.yStart < 0.0001) ||
		(axis % 3 == 2 && tree.zEnd - tree.zStart < 0.0001))
	{
		tree.leftChild = nullptr;
		tree.rightChild = nullptr;
		return;
	}

	AccelTreeNode *left = new AccelTreeNode(), *right = new AccelTreeNode();
	std::vector<int> *leftTri = new std::vector<int>(),
		*rightTri = new std::vector<int>(),
		*bothTri = new std::vector<int>();

	float mid;

	// calculate the best split out of 10 samples:
	mid = calcBestSplit(tree, axis);

	// split the space in 2 smaller subspaces and recursively call this function on those subspaces
	if (axis % 3 == 0) {
		//split x-axis

		// split the space in half
		(*left) = AccelTreeNode(tree.xStart, mid, tree.yStart, tree.yEnd, tree.zStart, tree.zEnd);
		(*right) = AccelTreeNode(mid, tree.xEnd, tree.yStart, tree.yEnd, tree.zStart, tree.zEnd);
	}
	else if (axis % 3 == 1) {
		//split y-axis

		// split the space in half
		(*left) = AccelTreeNode(tree.xStart, tree.xEnd, tree.yStart, mid, tree.zStart, tree.zEnd);
		(*right) = AccelTreeNode(tree.xStart, tree.xEnd, mid, tree.yEnd, tree.zStart, tree.zEnd);
	}
	else { //if (axis % 3 == 2) {
		   //split z-axis

		   // split the space in half
		(*left) = AccelTreeNode(tree.xStart, tree.xEnd, tree.yStart, tree.yEnd, tree.zStart, mid);
		(*right) = AccelTreeNode(tree.xStart, tree.xEnd, tree.yStart, tree.yEnd, mid, tree.zEnd);
	}

	// calculate which triangles should be in the left, right, or both parts.
	for (std::vector<int>::const_iterator it = tree.triangles.begin(); it < tree.triangles.end(); ++it)
	{
		Triangle t = MyMesh.triangles[(*it)];
		if (MyMesh.vertices[t.v[0]].p[axis % 3] < mid && MyMesh.vertices[t.v[1]].p[axis % 3] < mid &&
			MyMesh.vertices[t.v[2]].p[axis % 3] < mid)
			(*leftTri).push_back(*it);
		else if (MyMesh.vertices[t.v[0]].p[axis % 3] > mid && MyMesh.vertices[t.v[1]].p[axis % 3] > mid &&
			MyMesh.vertices[t.v[2]].p[axis % 3] > mid)
			(*rightTri).push_back(*it);
		else
			(*bothTri).push_back(*it);
	}

	// Assign the triangles to the right node
	tree.triangles = *bothTri;
	(*left).triangles = *leftTri;
	(*right).triangles = *rightTri;

	// Set the parent lists of the left/right child nodes
	(*left).parentList = tree.parentList;
	(*left).parentList.push_back(left);
	(*right).parentList = tree.parentList;
	(*right).parentList.push_back(right);

	// Set left/right nodes as child nodes
	tree.leftChild = left;
	tree.rightChild = right;

	// recursively calculate the KDtree subspaces of the two parts
	splitSpaces((*left), axis + 1);
	++treeNodes;
	splitSpaces((*right), axis + 1);
	++treeNodes;
}

float calcBestSplit(AccelTreeNode &tree, int axis)
{
	float mid;
	float result = 0;
	float cost = std::numeric_limits<float>::max(), tempCost;

	for (int i = 1; i < TREE_ACCURACY; ++i)
	{
		float leftSize, rightSize, totalSize;

		if (axis % 3 == 0) {
			//split x-axis
			mid = ((tree.xEnd - tree.xStart) * i / TREE_ACCURACY) + tree.xStart;
			leftSize = mid - tree.xStart;
			rightSize = tree.xEnd - mid;
			totalSize = tree.xEnd - tree.xStart;
		}
		else if (axis % 3 == 1) {
			//split y-axis
			mid = ((tree.yEnd - tree.yStart) * i / TREE_ACCURACY) + tree.yStart;
			leftSize = mid - tree.yStart;
			rightSize = tree.yEnd - mid;
			totalSize = tree.yEnd - tree.yStart;
		}
		else { //if (axis % 3 == 2) {
			   //split z-axis
			mid = ((tree.zEnd - tree.zStart) * i / TREE_ACCURACY) + tree.zStart;
			leftSize = mid - tree.zStart;
			rightSize = tree.zEnd - mid;
			totalSize = tree.zEnd - tree.zStart;
		}

		std::vector<int> leftTri, rightTri, bothTri;

		// calculate which triangles should be in the left, right, or both parts.
		for (std::vector<int>::const_iterator it = tree.triangles.begin(); it < tree.triangles.end(); ++it)
		{
			Triangle t = MyMesh.triangles[(*it)];
			if (MyMesh.vertices[t.v[0]].p[axis % 3] < mid && MyMesh.vertices[t.v[1]].p[axis % 3] < mid &&
				MyMesh.vertices[t.v[2]].p[axis % 3] < mid)
				leftTri.push_back(*it);
			else if (MyMesh.vertices[t.v[0]].p[axis % 3] > mid && MyMesh.vertices[t.v[1]].p[axis % 3] > mid &&
				MyMesh.vertices[t.v[2]].p[axis % 3] > mid)
				rightTri.push_back(*it);
			else
				bothTri.push_back(*it);
		}

		tempCost = bothTri.size() * totalSize * 20.0f / axis + leftTri.size() * leftSize + rightTri.size() * rightSize;

		if (tempCost < cost)
		{
			cost = tempCost;
			result = mid;
		}

	}

	return result;
}


// Recursively find in which child node of the given parent the given position is
// If you don't know a parent, use the treeRoot variable (it's the root of the tree) and use axis = 0
inline AccelTreeNode findChildNode(const AccelTreeNode &parent, int axis, const Vec3Df &position)
{
	// if this node doesn't have any children, we found the right node
	if (parent.leftChild == nullptr)
		return parent;

	// else continue looking
	if (axis % 3 == 0)
	{
		// check x-axis
		if ((*parent.leftChild).xEnd > position.p[0])
			return findChildNode((*parent.leftChild), axis + 1, position);
		else
			return findChildNode((*parent.rightChild), axis + 1, position);
	}
	else if (axis % 3 == 1)
	{
		// check y-axis
		if ((*parent.leftChild).yEnd > position.p[1])
			return findChildNode((*parent.leftChild), axis + 1, position);
		else
			return findChildNode((*parent.rightChild), axis + 1, position);
	}
	else
	{
		// check z-axis
		if ((*parent.leftChild).zEnd > position.p[2])
			return findChildNode((*parent.leftChild), axis + 1, position);
		else
			return findChildNode((*parent.rightChild), axis + 1, position);
	}
}

// Find the point where the ray hits the outside of the current node
inline AccelTreeNode findNextNode(const AccelTreeNode &curN, const Vec3Df &position, const Vec3Df &destination)
{

	Vec3Df dir = destination - position;
	float f, tempx, tempy, tempz;
	AccelTreeNode node = treeRoot;

	// Check the y-z faces of the space
	if (dir.p[0] != 0)
	{
		if (dir.p[0] < 0)
		{
			// negative direction:
			// check xStart
			f = (curN.xStart - position.p[0]) / dir.p[0];
			tempy = f * dir.p[1] + position.p[1];
			tempz = f * dir.p[2] + position.p[2];

			// check if succesful
			if (tempy > curN.yStart && tempy < curN.yEnd && tempz > curN.zStart && tempz < curN.zEnd)
			{
				// return the root if we would otherwise go outside of the root.
				if (curN.xStart == treeRoot.xStart)
					return treeRoot;

				// return the right node
				for (int axis = 0;; ++axis)
				{
					if (axis % 3 == 0)
					{
						// check x-axis
						if ((*node.leftChild).xEnd >= curN.xStart)
							node = *node.leftChild;
						else
							node = *node.rightChild;
					}
					else if (axis % 3 == 1)
					{
						// check y-axis
						if (tempy < (*node.leftChild).yEnd)
							node = *node.leftChild;
						else if (tempy > (*node.leftChild).yEnd || dir[1] > 0)
							node = *node.rightChild;
						else
							node = *node.leftChild;
					}
					else
					{
						// check z-axis
						if (tempz < (*node.leftChild).zEnd)
							node = *node.leftChild;
						else if (tempz > (*node.leftChild).zEnd || dir[2] > 0)
							node = *node.rightChild;
						else
							node = *node.leftChild;
					}
					if (node.leftChild == nullptr)
						return node;
				}
			}
		}
		else
		{
			// positive direction:
			// check xEnd
			f = (curN.xEnd - position.p[0]) / dir.p[0];
			tempy = f * dir.p[1] + position.p[1];
			tempz = f * dir.p[2] + position.p[2];

			// check if succesful
			if (tempy > curN.yStart && tempy < curN.yEnd && tempz > curN.zStart && tempz < curN.zEnd)
			{
				// return the root if we would otherwise go outside of the root.
				if (curN.xEnd == treeRoot.xEnd)
					return treeRoot;

				// return the right node
				for (int axis = 0;; ++axis)
				{
					if (axis % 3 == 0)
					{
						// check x-axis
						if ((*node.leftChild).xEnd <= curN.xEnd)
							node = *node.rightChild;
						else
							node = *node.leftChild;
					}
					else if (axis % 3 == 1)
					{
						// check y-axis
						if (tempy < (*node.leftChild).yEnd)
							node = *node.leftChild;
						else if (tempy > (*node.leftChild).yEnd || dir[1] > 0)
							node = *node.rightChild;
						else
							node = *node.leftChild;
					}
					else
					{
						// check z-axis
						if (tempz < (*node.leftChild).zEnd)
							node = *node.leftChild;
						else if (tempz > (*node.leftChild).zEnd || dir[2] > 0)
							node = *node.rightChild;
						else
							node = *node.leftChild;
					}
					if (node.leftChild == nullptr)
						return node;
				}
			}
		}
	}

	// Check the x-z faces of the space
	if (dir.p[1] != 0)
	{
		if (dir.p[1] < 0)
		{
			// negative direction:
			// check yStart
			f = (curN.yStart - position.p[1]) / dir.p[1];
			tempx = f * dir.p[0] + position.p[0];
			tempz = f * dir.p[2] + position.p[2];

			// check if succesful
			if (tempx > curN.xStart && tempx < curN.xEnd && tempz > curN.zStart && tempz < curN.zEnd)
			{
				// return the root if we would otherwise go outside of the root.
				if (curN.yStart == treeRoot.yStart)
					return treeRoot;

				// return the right node
				for (int axis = 0;; ++axis)
				{
					if (axis % 3 == 0)
					{
						// check x-axis
						if (tempx < (*node.leftChild).xEnd)
							node = *node.leftChild;
						else if (tempx > (*node.leftChild).xEnd || dir[0] > 0)
							node = *node.rightChild;
						else
							node = *node.leftChild;
					}
					else if (axis % 3 == 1)
					{
						// check y-axis
						if ((*node.leftChild).yEnd >= curN.yStart)
							node = *node.leftChild;
						else
							node = *node.rightChild;
					}
					else
					{
						// check z-axis
						if (tempz < (*node.leftChild).zEnd)
							node = *node.leftChild;
						else if (tempz > (*node.leftChild).zEnd || dir[2] > 0)
							node = *node.rightChild;
						else
							node = *node.leftChild;
					}
					if (node.leftChild == nullptr)
						return node;
				}
			}
		}
		else
		{
			// positive direction:
			// check yEnd
			f = (curN.yEnd - position.p[1]) / dir.p[1];
			tempx = f * dir.p[0] + position.p[0];
			tempz = f * dir.p[2] + position.p[2];

			// check if succesful
			if (tempx > curN.xStart && tempx < curN.xEnd && tempz > curN.zStart && tempz < curN.zEnd)
			{
				// return the root if we would otherwise go outside of the root.
				if (curN.yEnd == treeRoot.yEnd)
					return treeRoot;

				// return the right node
				for (int axis = 0;; ++axis)
				{
					if (axis % 3 == 0)
					{
						// check x-axis
						if (tempx < (*node.leftChild).xEnd)
							node = *node.leftChild;
						else if (tempx > (*node.leftChild).xEnd || dir[0] > 0)
							node = *node.rightChild;
						else
							node = *node.leftChild;
					}
					else if (axis % 3 == 1)
					{
						// check y-axis
						if ((*node.leftChild).yEnd <= curN.yEnd)
							node = *node.rightChild;
						else
							node = *node.leftChild;
					}
					else
					{
						// check z-axis
						if (tempz < (*node.leftChild).zEnd)
							node = *node.leftChild;
						else if (tempz > (*node.leftChild).zEnd || dir[2] > 0)
							node = *node.rightChild;
						else
							node = *node.leftChild;
					}
					if (node.leftChild == nullptr)
						return node;
				}
			}
		}
	}

	// Check the x-y faces of the space
	if (dir.p[2] != 0)
	{
		if (dir.p[2] < 0)
		{
			// negative direction:
			// check yStart
			f = (curN.zStart - position.p[2]) / dir.p[2];
			tempx = f * dir.p[0] + position.p[0];
			tempy = f * dir.p[1] + position.p[1];

			// return the root if we would otherwise go outside of the root.
			if (curN.zStart == treeRoot.zStart)
				return treeRoot;

			// return the right node
			for (int axis = 0;; ++axis)
			{
				if (axis % 3 == 0)
				{
					// check x-axis
					if (tempx < (*node.leftChild).xEnd)
						node = *node.leftChild;
					else if(tempx > (*node.leftChild).xEnd || dir[0] > 0)
						node = *node.rightChild;
					else
						node = *node.leftChild;
				}
				else if (axis % 3 == 1)
				{
					// check y-axis
					if (tempy < (*node.leftChild).yEnd)
						node = *node.leftChild;
					else if (tempy > (*node.leftChild).yEnd || dir[1] > 0)
						node = *node.rightChild;
					else
						node = *node.leftChild;
				}
				else
				{
					// check z-axis
					if ((*node.leftChild).zEnd >= curN.zStart)
						node = *node.leftChild;
					else
						node = *node.rightChild;
				}
				if (node.leftChild == nullptr)
					return node;
			}
		}
		else
		{
			// positive direction:
			// check yEnd
			f = (curN.zEnd - position.p[2]) / dir.p[2];
			tempx = f * dir.p[0] + position.p[0];
			tempy = f * dir.p[1] + position.p[1];

			// return the root if we would otherwise go outside of the root.
			if (curN.zEnd == treeRoot.zEnd)
				return treeRoot;

			// return the right node
			for (int axis = 0;; ++axis)
			{
				if (axis % 3 == 0)
				{
					// check x-axis
					if (tempx < (*node.leftChild).xEnd)
						node = *node.leftChild;
					else if (tempx > (*node.leftChild).xEnd || dir[0] > 0)
						node = *node.rightChild;
					else
						node = *node.leftChild;
				}
				else if (axis % 3 == 1)
				{
					// check y-axis
					if (tempy < (*node.leftChild).yEnd)
						node = *node.leftChild;
					else if (tempy > (*node.leftChild).yEnd || dir[1] > 0)
						node = *node.rightChild;
					else
						node = *node.leftChild;
				}
				else
				{
					// check z-axis
					if ((*node.leftChild).zEnd <= curN.zEnd)
						node = *node.rightChild;
					else
						node = *node.leftChild;
				}
				if (node.leftChild == nullptr)
					return node;
			}
		}
	}

	return treeRoot;
}

inline bool contains(const std::vector<AccelTreeNode*> &vec, const AccelTreeNode &element)
{
	for (std::vector<AccelTreeNode*>::const_iterator it = vec.begin(); it != vec.end(); ++it)
		if ((**it).xStart == element.xStart &&
			(**it).yStart == element.yStart &&
			(**it).zStart == element.zStart &&
			(**it).xEnd == element.xEnd &&
			(**it).yEnd == element.yEnd &&
			(**it).zEnd == element.zEnd)
		{
			return true;
		}
	return false;
}

inline void projectOriginOnRoot(Vec3Df &origin, Vec3Df &dest)
{
	Vec3Df hitRoot = calculateProjectionOnRoot(origin, dest);

	// avoid floating point error
	// !!! NEEDS BETTER SOLUTION !!!
	hitRoot += (dest / 1000.f);

	dest += hitRoot - origin;
	origin = hitRoot;
}

// If the camera is outside of the root, use this to move the origin inside of the root.
// this function looks a lot like the findNodeBoxHitPoint function
inline Vec3Df calculateProjectionOnRoot(Vec3Df &origin, Vec3Df &dest)
{
	float i, tempx, tempy, tempz;

	// check y-z faces of the root
	if (origin.p[0] != dest.p[0] && !(origin.p[0] > treeRoot.xStart && origin.p[0] < treeRoot.xEnd))
	{
		if (std::abs(treeRoot.xStart - origin.p[0]) < std::abs(treeRoot.xEnd - origin.p[0]))
		{
			// closer to xStart
			tempx = treeRoot.xStart - origin.p[0];
			i = tempx / (dest.p[0] - origin.p[0]);

			tempy = i * (dest.p[1] - origin.p[1]) + origin.p[1];
			tempz = i * (dest.p[2] - origin.p[2]) + origin.p[2];

			// check if succesful
			if (tempy > treeRoot.yStart && tempy < treeRoot.yEnd && tempz > treeRoot.zStart && tempz < treeRoot.zEnd)
				return Vec3Df(treeRoot.xStart, tempy, tempz);
		}
		else
		{
			// closer to xEnd
			tempx = treeRoot.xEnd - origin.p[0];
			i = tempx / (dest.p[0] - origin.p[0]);

			tempy = i * (dest.p[1] - origin.p[1]) + origin.p[1];
			tempz = i * (dest.p[2] - origin.p[2]) + origin.p[2];

			// check if succesful
			if (tempy > treeRoot.yStart && tempy < treeRoot.yEnd && tempz > treeRoot.zStart && tempz < treeRoot.zEnd)
				return Vec3Df(treeRoot.xEnd, tempy, tempz);
		}
	}

	// check x-z faces of the root
	if (origin.p[1] != dest.p[1] && !(origin.p[1] > treeRoot.yStart && origin.p[1] < treeRoot.yEnd))
	{
		if (std::abs(treeRoot.yStart - origin.p[1]) < std::abs(treeRoot.yEnd - origin.p[1]))
		{
			// closer to yStart
			tempy = treeRoot.yStart - origin.p[1];
			i = tempy / (dest.p[1] - origin.p[1]);

			tempx = i * (dest.p[0] - origin.p[0]) + origin.p[0];
			tempz = i * (dest.p[2] - origin.p[2]) + origin.p[2];

			// check if succesful
			if (tempx > treeRoot.xStart && tempx < treeRoot.xEnd && tempz > treeRoot.zStart && tempz < treeRoot.zEnd)
				return Vec3Df(tempx, treeRoot.yStart, tempz);
		}
		else
		{
			// closer to yEnd
			tempy = treeRoot.yEnd - origin.p[1];
			i = tempy / (dest.p[1] - origin.p[1]);

			tempx = i * (dest.p[0] - origin.p[0]) + origin.p[0];
			tempz = i * (dest.p[2] - origin.p[2]) + origin.p[2];

			// check if succesful
			if (tempx > treeRoot.xStart && tempx < treeRoot.xEnd && tempz > treeRoot.zStart && tempz < treeRoot.zEnd)
				return Vec3Df(tempx, treeRoot.yEnd, tempz);
		}
	}

	// check x-y faces of the root
	if (origin.p[2] != dest.p[2] && !(origin.p[2] > treeRoot.zStart && origin.p[2] < treeRoot.zEnd))
	{
		if (std::abs(treeRoot.zStart - origin.p[2]) < std::abs(treeRoot.zEnd - origin.p[2]))
		{
			// closer to zStart
			tempz = treeRoot.zStart - origin.p[2];
			i = tempz / (dest.p[2] - origin.p[2]);

			if (tempx > treeRoot.xStart && tempx < treeRoot.xEnd && tempy > treeRoot.yStart && tempy < treeRoot.yEnd)
				return Vec3Df(tempx, tempy, treeRoot.zStart);
		}
		else
		{
			// closer to zEnd
			tempz = treeRoot.zEnd - origin.p[2];
			i = tempz / (dest.p[2] - origin.p[2]);

			tempx = i * (dest.p[0] - origin.p[0]) + origin.p[0];
			tempy = i * (dest.p[1] - origin.p[1]) + origin.p[1];

			if (tempx > treeRoot.xStart && tempx < treeRoot.xEnd && tempy > treeRoot.yStart && tempy < treeRoot.yEnd)
				return Vec3Df(tempx, tempy, treeRoot.zEnd);
		}
	}

	return origin;
}