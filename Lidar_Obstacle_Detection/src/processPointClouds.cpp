// PCL lib Functions for processing point clouds 

#include <unordered_set>
#include "processPointClouds.h"


//constructor:
template<typename PointT>
ProcessPointClouds<PointT>::ProcessPointClouds() {}


//de-constructor:
template<typename PointT>
ProcessPointClouds<PointT>::~ProcessPointClouds() {}


template<typename PointT>
void ProcessPointClouds<PointT>::numPoints(typename pcl::PointCloud<PointT>::Ptr cloud)
{
    std::cout << cloud->points.size() << std::endl;
}


// Filter the input point cloud.
// -- Down samples the point cloud based on the resolution filterRes.
// -- Crops the point cloud to used the points within interested areas.
template<typename PointT>
typename pcl::PointCloud<PointT>::Ptr ProcessPointClouds<PointT>::FilterCloud(
    typename pcl::PointCloud<PointT>::Ptr cloud, 
    float filterRes, 
    Eigen::Vector4f minPoint, 
    Eigen::Vector4f maxPoint)
{

    // Time segmentation process
    auto startTime = std::chrono::steady_clock::now();

    // Down sample the point cloud
    pcl::VoxelGrid<PointT> vg;
    vg.setInputCloud (cloud);
    vg.setLeafSize (filterRes, filterRes, filterRes);
    vg.filter (*cloud);

    // Create a crop box
    pcl::CropBox<PointT> region(true);
    region.setInputCloud(cloud);
    region.setMin(minPoint);
    region.setMax(maxPoint);
    region.filter(*cloud);

    // Filter ego car roof voxels
    std::vector<int> indices; 
    pcl::CropBox<PointT> roof(true);
    roof.setInputCloud(cloud);
    roof.setMin(Eigen::Vector4f(-1.5, -1.7, -1, 1));
    roof.setMax(Eigen::Vector4f(2.6, 1.7, -.4, 1));
    roof.filter(indices);

    pcl::PointIndices::Ptr inliers {new pcl::PointIndices};
    for (int index : indices) 
        inliers->indices.push_back(index);
    
    pcl::ExtractIndices<PointT> ec;
    ec.setInputCloud(cloud);
    ec.setIndices(inliers);
    ec.setNegative(true);
    ec.filter(*cloud);
    
    auto endTime = std::chrono::steady_clock::now();
    auto elapsedTime = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);
    std::cout << "filtering took " << elapsedTime.count() << " milliseconds" << std::endl;

    return cloud;

}

// Segment street point cloud and obstacle point cloud from the input point cloud.
template<typename PointT>
std::pair<typename pcl::PointCloud<PointT>::Ptr, typename pcl::PointCloud<PointT>::Ptr> 
    ProcessPointClouds<PointT>::SeparateClouds(
        pcl::PointIndices::Ptr inliers, 
        typename pcl::PointCloud<PointT>::Ptr cloud) 
{
    typename pcl::PointCloud<PointT>::Ptr cloud_road (new pcl::PointCloud<PointT>);
    typename pcl::PointCloud<PointT>::Ptr cloud_obstacle (new pcl::PointCloud<PointT>);

    for (int index : inliers->indices) 
        cloud_road->points.push_back(cloud->points[index]);

    pcl::ExtractIndices<PointT> extract;
    extract.setInputCloud(cloud);
    extract.setIndices(inliers);
    extract.setNegative(true);
    extract.filter(*cloud_obstacle);

    std::pair<typename pcl::PointCloud<PointT>::Ptr, typename pcl::PointCloud<PointT>::Ptr> segResult(cloud_obstacle, cloud_road);
    return segResult;
}

// 3D RANSAC implementation
template<typename PointT>
pcl::PointIndices::Ptr ProcessPointClouds<PointT>::Ransac(
    typename pcl::PointCloud<PointT>::Ptr cloud, 
    int maxIterations, 
    float distanceTol)
{
	pcl::PointIndices::Ptr inliersResult {new pcl::PointIndices};
	srand(time(NULL));
	
	for (int i = 0; i < maxIterations; i++){
        // Randomly pick 3 unique points
        std::unordered_set<int> inliers;
        while (inliers.size() < 3)
            inliers.insert(rand() % (cloud->points.size()));
        
        // Generate coefficients for the plane that contains these 3 points.
        auto iter = inliers.begin();

		double x1 = cloud->points[*iter].x;
		double y1 = cloud->points[*iter].y;
        double z1 = cloud->points[*iter].z;
		iter++;
		double x2 = cloud->points[*iter].x;
		double y2 = cloud->points[*iter].y;
		double z2 = cloud->points[*iter].z;
        iter++;
		double x3 = cloud->points[*iter].x;
		double y3 = cloud->points[*iter].y;
		double z3 = cloud->points[*iter].z;

		double a = (y2-y1)*(z3-z1)-(z2-z1)*(y3-y1);
		double b = (z2-z1)*(x3-x1)-(x2-x1)*(z3-z1);
		double c = (x1-x1)*(y3-y1)-(y2-y1)*(x3-x1);
		double d = -(a*x1+b*y1+c*z1);

		pcl::PointIndices::Ptr currentInliers {new pcl::PointIndices};
		for (int index = 0; index < cloud->points.size(); index++) {
			if ( fabs(a * cloud->points[index].x + b * cloud->points[index].y + c * cloud->points[index].z + d) / sqrt(a * a + b * b + c * c) <= distanceTol) {
				currentInliers->indices.push_back(index);
			}
		}
		if (currentInliers->indices.size() > inliersResult->indices.size()){
			inliersResult = currentInliers;
		}

	}
	
	return inliersResult;

}

// Segment plane from the input point cloud
template<typename PointT>
std::pair<typename pcl::PointCloud<PointT>::Ptr, typename pcl::PointCloud<PointT>::Ptr> 
    ProcessPointClouds<PointT>::SegmentPlane(typename pcl::PointCloud<PointT>::Ptr cloud, int maxIterations, float distanceThreshold)
{
    // Time segmentation process
    auto startTime = std::chrono::steady_clock::now();

    pcl::PointIndices::Ptr inliers = Ransac(cloud, maxIterations, distanceThreshold);

    /* 
     * Below are the original pcl based SegmentPlane.
     * // Create the segmentation object
     * pcl::SACSegmentation<PointT> seg;
     * pcl::PointIndices::Ptr inliers {new pcl::PointIndices};
     * pcl::ModelCoefficients::Ptr coefficients {new pcl::ModelCoefficients};
     *
     * // Optional
     * seg.setOptimizeCoefficients(true);
     * // Mandatory
     * seg.setModelType(pcl::SACMODEL_PLANE);
     * seg.setMethodType(pcl::SAC_RANSAC);
     * seg.setMaxIterations(maxIterations);
     * seg.setDistanceThreshold(distanceThreshold);
     * 
     * seg.setInputCloud(cloud);
     * seg.segment(*inliers, *coefficients);
    */
    
    if (inliers->indices.size() == 0) {
        std::cout << "Could not estimate a planar model for the given dataset." << std::endl;
    }
    

    auto endTime = std::chrono::steady_clock::now();
    auto elapsedTime = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);
    std::cout << "plane segmentation took " << elapsedTime.count() << " milliseconds" << std::endl;

    std::pair<typename pcl::PointCloud<PointT>::Ptr, typename pcl::PointCloud<PointT>::Ptr> segResult = SeparateClouds(inliers, cloud);
    return segResult;
}

template<typename PointT>
void ProcessPointClouds<PointT>::proximity(
    typename pcl::PointCloud<PointT>::Ptr points, 
    int index, 
    pcl::PointIndices* cluster_indices, 
    std::set<int>* processed, 
    KdTree* tree, 
    float distanceTol) 
{
	if (processed->find(index) == processed->end())
	{
		processed->insert(index);
		cluster_indices->indices.push_back(index);
		std::vector<int> nearby_point;
		nearby_point = tree->search({(*points)[index].x, (*points)[index].y, (*points)[index].z}, distanceTol);
		for (int index : nearby_point)
		{
			proximity(points, index, cluster_indices, processed, tree, distanceTol);
		}
	}

}

// Euclidean clustering implementation.
template<typename PointT>
std::vector<pcl::PointIndices> ProcessPointClouds<PointT>::euclideanCluster(
    typename pcl::PointCloud<PointT>::Ptr points, 
    KdTree* tree, 
    float distanceTol, 
    int minClusterSize, 
    int maxClusterSize)
{
	std::vector<pcl::PointIndices> clusters;
	std::set<int> processed;

	for (int index = 0; index < points->size(); index++) 
	{
		if (processed.find(index) == processed.end())
		{
			pcl::PointIndices cluster;
			proximity(points, index, &cluster, &processed, tree, distanceTol);
            if (cluster.indices.size() >= minClusterSize && cluster.indices.size() <= maxClusterSize)
			    clusters.push_back(cluster);
		}
	}
	
	return clusters;

}

template<typename PointT>
std::vector<typename pcl::PointCloud<PointT>::Ptr> ProcessPointClouds<PointT>::Clustering(
    typename pcl::PointCloud<PointT>::Ptr cloud, 
    float clusterTolerance, 
    int minSize, 
    int maxSize)
{

    // Time clustering process
    auto startTime = std::chrono::steady_clock::now();

    std::vector<typename pcl::PointCloud<PointT>::Ptr> clusters;

    // Clustering use self-implemented euclidean cluster method.
    KdTree* tree  {new KdTree};
    for (int i = 0; i < cloud->size(); i++)
        tree->insert({(*cloud)[i].x, (*cloud)[i].y, (*cloud)[i].z}, i);

    std::vector<pcl::PointIndices> cluster_indices = euclideanCluster(cloud, tree, clusterTolerance, minSize, maxSize);

    /*
     * Below are the original pcl based Clustering.
     * typename pcl::search::KdTree<PointT>::Ptr tree (new pcl::search::KdTree<PointT>);
     * tree->setInputCloud (cloud);
     * 
     * std::vector<pcl::PointIndices> cluster_indices;
     * pcl::EuclideanClusterExtraction<PointT> ec;
     * ec.setClusterTolerance(clusterTolerance); 
     * ec.setMinClusterSize(minSize);
     * ec.setMaxClusterSize(maxSize);
     * ec.setSearchMethod(tree);
     * ec.setInputCloud(cloud);
     * ec.extract(cluster_indices);
    */
    
    for (pcl::PointIndices getIndices : cluster_indices)
    {
        typename pcl::PointCloud<PointT>::Ptr cluster (new pcl::PointCloud<PointT>);
        for (int index : getIndices.indices)
        {
            cluster->points.push_back(cloud->points[index]);
        }
        cluster->width = cluster->points.size();
        cluster->height = 1;
        cluster->is_dense = true;

        clusters.push_back(cluster);
    }
  

    auto endTime = std::chrono::steady_clock::now();
    auto elapsedTime = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);
    std::cout << "clustering took " << elapsedTime.count() << " milliseconds and found " << clusters.size() << " clusters" << std::endl;

    return clusters;
}


template<typename PointT>
Box ProcessPointClouds<PointT>::BoundingBox(typename pcl::PointCloud<PointT>::Ptr cluster)
{

    // Find bounding box for one of the clusters
    PointT minPoint, maxPoint;
    pcl::getMinMax3D(*cluster, minPoint, maxPoint);

    Box box;
    box.x_min = minPoint.x;
    box.y_min = minPoint.y;
    box.z_min = minPoint.z;
    box.x_max = maxPoint.x;
    box.y_max = maxPoint.y;
    box.z_max = maxPoint.z;

    return box;
}

// Create bounding box based on PCA.
template<typename PointT>
BoxQ ProcessPointClouds<PointT>::BoundingBoxQ(typename pcl::PointCloud<PointT>::Ptr cluster)
{

    // Find bounding box for one of the clusters
    // Compute principal directions
    Eigen::Vector4f pcaCentroid;
    pcl::compute3DCentroid(*cluster, pcaCentroid);
    Eigen::Matrix3f covariance;
    computeCovarianceMatrixNormalized(*cluster, pcaCentroid, covariance);
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3f> eigen_solver(covariance, Eigen::ComputeEigenvectors);
    Eigen::Matrix3f eigenVectorsPCA = eigen_solver.eigenvectors();

    // This line is necessary for proper orientation in some cases. The numbers come out the same without it, but
    // the signs are different and the box doesn't get correctly oriented in some cases.
    eigenVectorsPCA.col(2) = eigenVectorsPCA.col(0).cross(eigenVectorsPCA.col(1));  

    Eigen::Matrix4f projectionTransform(Eigen::Matrix4f::Identity());
    projectionTransform.block<3,3>(0,0) = eigenVectorsPCA.transpose();
    projectionTransform.block<3,1>(0,3) = -1.f * (projectionTransform.block<3,3>(0,0) * pcaCentroid.head<3>());
    typename pcl::PointCloud<PointT>::Ptr cloudPointsProjected (new pcl::PointCloud<PointT>);
    pcl::transformPointCloud(*cluster, *cloudPointsProjected, projectionTransform);
    
    // Get the minimum and maximum points of the transformed cloud.
    PointT minPoint, maxPoint;
    pcl::getMinMax3D(*cloudPointsProjected, minPoint, maxPoint);
    const Eigen::Vector3f meanDiagonal = 0.5f*(maxPoint.getVector3fMap() + minPoint.getVector3fMap());
    
    // Final transform
    const Eigen::Quaternionf bboxQuaternion(eigenVectorsPCA); 
    const Eigen::Vector3f bboxTransform = eigenVectorsPCA * meanDiagonal + pcaCentroid.head<3>();

    BoxQ box;
    box.bboxTransform = bboxTransform;
    box.bboxQuaternion = bboxQuaternion;
    box.cube_length = maxPoint.x - minPoint.x;
    box.cube_width = maxPoint.y - minPoint.y;
    box.cube_height = maxPoint.z - minPoint.z;

    return box;

}

template<typename PointT>
void ProcessPointClouds<PointT>::savePcd(typename pcl::PointCloud<PointT>::Ptr cloud, std::string file)
{
    pcl::io::savePCDFileASCII (file, *cloud);
    std::cerr << "Saved " << cloud->points.size () << " data points to "+file << std::endl;
}


template<typename PointT>
typename pcl::PointCloud<PointT>::Ptr ProcessPointClouds<PointT>::loadPcd(std::string file)
{

    typename pcl::PointCloud<PointT>::Ptr cloud (new pcl::PointCloud<PointT>);

    if (pcl::io::loadPCDFile<PointT> (file, *cloud) == -1) //* load the file
    {
        PCL_ERROR ("Couldn't read file \n");
    }
    std::cerr << "Loaded " << cloud->points.size () << " data points from "+file << std::endl;

    return cloud;
}


template<typename PointT>
std::vector<boost::filesystem::path> ProcessPointClouds<PointT>::streamPcd(std::string dataPath)
{

    std::vector<boost::filesystem::path> paths(boost::filesystem::directory_iterator{dataPath}, boost::filesystem::directory_iterator{});

    // sort files in accending order so playback is chronological
    sort(paths.begin(), paths.end());

    return paths;

}