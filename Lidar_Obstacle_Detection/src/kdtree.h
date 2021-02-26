// Author: Junjie Yan
// KD Tree implementation

#include <iostream>
#include <vector>
#include <string>


// Represent the node of a kd tree
struct Node
{
	std::vector<float> point;
	int id;
	Node* left;
	Node* right;

	Node(std::vector<float> arr, int setId)
	:	point(arr), id(setId), left(NULL), right(NULL)
	{}
};

struct KdTree
{
	Node* root;

	KdTree()
	: root(NULL)
	{}

	void insertHelper(Node **node, std::vector<float> point, int id, int depth)
	{
		if (*node == NULL)
		{
			*node = new Node(point, id);
		}
		else 
		{
			// caculate current dim
			uint cd = depth % point.size();
			if (point[cd] < (*node)->point[cd])
				insertHelper(&((*node)->left), point, id, depth + 1);
			else
				insertHelper(&((*node)->right), point, id, depth + 1);
		}
	}

	void insert(std::vector<float> point, int id)
	{
		insertHelper(&root, point, id, 0);
	}

	void searchHelper(Node* node, std::vector<float> target, float distanceTol, std::vector<int>& ids, int depth)
	{
		if (node == NULL)
			return;
		if ((node->point[0]>= target[0] - distanceTol) && 
			(node->point[0]<= target[0] + distanceTol) && 
			(node->point[1]>= target[1] - distanceTol) && 
			(node->point[1]<= target[1] + distanceTol) &&
			(node->point[2]>= target[2] - distanceTol) && 
			(node->point[2]<= target[2] + distanceTol))
			{
				ids.push_back(node->id);
			}
		if (target[depth % target.size()] - distanceTol < node->point[depth % target.size()])
		{
			searchHelper(node->left, target, distanceTol, ids, depth + 1);
		}
		if (target[depth % target.size()] + distanceTol > node->point[depth % target.size()])
		{
			searchHelper(node->right, target, distanceTol, ids, depth + 1);
		}
	}

	// return a list of point ids in the tree that are within distance of target
	std::vector<int> search(std::vector<float> target, float distanceTol)
	{
		std::vector<int> ids;
		searchHelper(root, target, distanceTol, ids, 0);
		return ids;
	}
	

};




