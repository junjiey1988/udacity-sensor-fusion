// Author: Junjie Yan
// Lidar Obstacle Detection Project
// Use self implemented RANSAC, KD-Tree and Euclidean clustering algorithem.

#include "render/render.h"
#include "processPointClouds.h"
#include "processPointClouds.cpp"

// Open 3D viewer and display the City Block given the PCD date.
void cityBlock(
    pcl::visualization::PCLVisualizer::Ptr& viewer, 
    ProcessPointClouds<pcl::PointXYZI>* pointProcessorI, 
    const pcl::PointCloud<pcl::PointXYZI>::Ptr& inputCloud)
{
    pcl::PointCloud<pcl::PointXYZI>::Ptr filterCloud = 
        pointProcessorI->FilterCloud(
            inputCloud, 
            0.1, 
            Eigen::Vector4f(-10, -6, -2, 1), 
            Eigen::Vector4f( 30, 6, 2, 1));

    // Segment road and obstacle
    std::pair<
        pcl::PointCloud<pcl::PointXYZI>::Ptr, 
        pcl::PointCloud<pcl::PointXYZI>::Ptr> segmentCloud = 
            pointProcessorI->SegmentPlane(inputCloud, 50, 0.2);

    // Render the street point cloud.
    renderPointCloud(viewer,segmentCloud.second,"planeCloud",Color(0,1,0));

    // Cluster the obstacle point cloud.
    std::vector<pcl::PointCloud<pcl::PointXYZI>::Ptr> cloudClusters = 
        pointProcessorI->Clustering(segmentCloud.first, 0.3, 50, 10000);

    int clusterId = 0;
    std::vector<Color> colors = {Color(1,0,0), Color(1,1,0), Color(0,0,1)};

    for(pcl::PointCloud<pcl::PointXYZI>::Ptr cluster : cloudClusters)
    {
        std::cout << "cluster size ";
        pointProcessorI->numPoints(cluster);
        renderPointCloud(
            viewer,cluster,
            "obstCloud"+std::to_string(clusterId),
            colors[clusterId % 3]);

        Box box = pointProcessorI->BoundingBox(cluster);
        // BoxQ is the bounding box generated after using PCA.
        // BoxQ box = pointProcessorI->BoundingBoxQ(cluster);
        renderBox(viewer, box, clusterId);

        ++clusterId;
    }

    


}


//setAngle: SWITCH CAMERA ANGLE {XY, TopDown, Side, FPS}
void initCamera(CameraAngle setAngle, pcl::visualization::PCLVisualizer::Ptr& viewer)
{

    viewer->setBackgroundColor (0, 0, 0);
    
    // set camera position and angle
    viewer->initCameraParameters();
    // distance away in meters
    int distance = 16;
    
    switch(setAngle)
    {
        case XY : 
            viewer->setCameraPosition(-distance, -distance, distance, 1, 1, 0); 
            break;
        case TopDown : 
            viewer->setCameraPosition(0, 0, distance, 1, 0, 1); 
            break;
        case Side : 
            viewer->setCameraPosition(0, -distance, 0, 0, 0, 1); 
            break;
        case FPS : 
            viewer->setCameraPosition(-10, 0, 0, 0, 0, 1);
    }

    if(setAngle!=FPS)
        viewer->addCoordinateSystem (1.0);
}


int main (int argc, char** argv)
{
    std::cout << "starting enviroment" << std::endl;

    pcl::visualization::PCLVisualizer::Ptr viewer (new pcl::visualization::PCLVisualizer ("3D Viewer"));
    CameraAngle setAngle = FPS;
    initCamera(setAngle, viewer);

    ProcessPointClouds<pcl::PointXYZI>* pointProcessorI (new ProcessPointClouds<pcl::PointXYZI>());
    std::vector<boost::filesystem::path> stream = pointProcessorI->streamPcd("../src/sensors/data/pcd/data_2");
    auto streamIterator = stream.begin();
    pcl::PointCloud<pcl::PointXYZI>::Ptr inputCloudI;
    
    while (!viewer->wasStopped ())
    {
        // Clear viewer
        viewer->removeAllPointClouds();
        viewer->removeAllShapes();

        // Load pcd and run obstacle detection process
        inputCloudI = pointProcessorI->loadPcd((*streamIterator).string());
        cityBlock(viewer, pointProcessorI, inputCloudI);

        streamIterator++;
        if(streamIterator == stream.end())
            streamIterator = stream.begin();

        viewer->spinOnce ();
    } 
}