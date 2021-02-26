
#include <iostream>
#include <algorithm>
#include <numeric>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

#include "camFusion.hpp"
#include "dataStructures.h"

using namespace std;


// Create groups of Lidar points whose projection into the camera falls into the same bounding box
void clusterLidarWithROI(std::vector<BoundingBox> &boundingBoxes, std::vector<LidarPoint> &lidarPoints, float shrinkFactor, cv::Mat &P_rect_xx, cv::Mat &R_rect_xx, cv::Mat &RT)
{
    // loop over all Lidar points and associate them to a 2D bounding box
    cv::Mat X(4, 1, cv::DataType<double>::type);
    cv::Mat Y(3, 1, cv::DataType<double>::type);

    for (auto it1 = lidarPoints.begin(); it1 != lidarPoints.end(); ++it1)
    {
        // assemble vector for matrix-vector-multiplication
        X.at<double>(0, 0) = it1->x;
        X.at<double>(1, 0) = it1->y;
        X.at<double>(2, 0) = it1->z;
        X.at<double>(3, 0) = 1;

        // project Lidar point into camera
        Y = P_rect_xx * R_rect_xx * RT * X;
        cv::Point pt;
        // pixel coordinates
        pt.x = Y.at<double>(0, 0) / Y.at<double>(2, 0); 
        pt.y = Y.at<double>(1, 0) / Y.at<double>(2, 0); 

        vector<vector<BoundingBox>::iterator> enclosingBoxes; // pointers to all bounding boxes which enclose the current Lidar point
        for (vector<BoundingBox>::iterator it2 = boundingBoxes.begin(); it2 != boundingBoxes.end(); ++it2)
        {
            // shrink current bounding box slightly to avoid having too many outlier points around the edges
            cv::Rect smallerBox;
            smallerBox.x = (*it2).roi.x + shrinkFactor * (*it2).roi.width / 2.0;
            smallerBox.y = (*it2).roi.y + shrinkFactor * (*it2).roi.height / 2.0;
            smallerBox.width = (*it2).roi.width * (1 - shrinkFactor);
            smallerBox.height = (*it2).roi.height * (1 - shrinkFactor);

            // check wether point is within current bounding box
            if (smallerBox.contains(pt))
            {
                enclosingBoxes.push_back(it2);
            }

        } // eof loop over all bounding boxes

        // check wether point has been enclosed by one or by multiple boxes
        if (enclosingBoxes.size() == 1)
        { 
            // add Lidar point to bounding box
            enclosingBoxes[0]->lidarPoints.push_back(*it1);
        }

    } // eof loop over all Lidar points
}

bool lidarCompare(const LidarPoint &lhs, const LidarPoint &rhs)  
{
    return lhs.x < rhs.x;
}

/* 
* The show3DObjects() function below can handle different output image sizes, but the text output has been manually tuned to fit the 2000x2000 size. 
* However, you can make this function work for other sizes too.
* For instance, to use a 1000x1000 size, adjusting the text positions by dividing them by 2.
*/
void show3DObjects(std::vector<BoundingBox> &boundingBoxes, cv::Size worldSize, cv::Size imageSize, bool bWait)
{
    // create topview image
    cv::Mat topviewImg(imageSize, CV_8UC3, cv::Scalar(255, 255, 255));

    for(auto it1=boundingBoxes.begin(); it1!=boundingBoxes.end(); ++it1)
    {
        // create randomized color for current 3D object
        cv::RNG rng(it1->boxID);
        cv::Scalar currColor = cv::Scalar(rng.uniform(0,150), rng.uniform(0, 150), rng.uniform(0, 150));

        // plot Lidar points into top view image
        int top=1e8, left=1e8, bottom=0.0, right=0.0; 
        float xwmin=1e8, ywmin=1e8, ywmax=-1e8;
        for (auto it2 = it1->lidarPoints.begin(); it2 != it1->lidarPoints.end(); ++it2)
        {
            // world coordinates
            float xw = (*it2).x; // world position in m with x facing forward from sensor
            float yw = (*it2).y; // world position in m with y facing left from sensor
            xwmin = xwmin<xw ? xwmin : xw;
            ywmin = ywmin<yw ? ywmin : yw;
            ywmax = ywmax>yw ? ywmax : yw;

            // top-view coordinates
            int y = (-xw * imageSize.height / worldSize.height) + imageSize.height;
            int x = (-yw * imageSize.width / worldSize.width) + imageSize.width / 2;

            // find enclosing rectangle
            top = top<y ? top : y;
            left = left<x ? left : x;
            bottom = bottom>y ? bottom : y;
            right = right>x ? right : x;

            // draw individual point
            cv::circle(topviewImg, cv::Point(x, y), 4, currColor, -1);
        }

        // draw enclosing rectangle
        cv::rectangle(topviewImg, cv::Point(left, top), cv::Point(right, bottom),cv::Scalar(0,0,0), 2);

        // augment object with some key data
        char str1[200], str2[200];
        sprintf(str1, "id=%d, #pts=%d", it1->boxID, (int)it1->lidarPoints.size());
        putText(topviewImg, str1, cv::Point2f(left-250, bottom+50), cv::FONT_ITALIC, 2, currColor);
        sprintf(str2, "xmin=%2.2f m, yw=%2.2f m", xwmin, ywmax-ywmin);
        putText(topviewImg, str2, cv::Point2f(left-250, bottom+125), cv::FONT_ITALIC, 2, currColor);  
        
        // Annotate 5 percentile point
        if (it1->lidarPoints.size() > 0)
        {
            auto five_per_it = it1->lidarPoints.begin() + it1->lidarPoints.size() / 20;
            std::nth_element(it1->lidarPoints.begin(), five_per_it, it1->lidarPoints.end(), lidarCompare);
            int y = (-five_per_it->x * imageSize.height / worldSize.height) + imageSize.height;
            int x = (-five_per_it->y * imageSize.width / worldSize.width) + imageSize.width / 2;
            cv::circle(topviewImg, cv::Point(x, y), 8, cv::Scalar(0, 0, 255), -1);
            char str3[200];
            sprintf(str2, "5 percentile x=%2.2f m", five_per_it->x);
            putText(topviewImg, str2, cv::Point2f(left-250, bottom+200), cv::FONT_ITALIC, 2, currColor);  
        }
    }

    // plot distance markers
    float lineSpacing = 2.0; // gap between distance markers
    int nMarkers = floor(worldSize.height / lineSpacing);
    for (size_t i = 0; i < nMarkers; ++i)
    {
        int y = (-(i * lineSpacing) * imageSize.height / worldSize.height) + imageSize.height;
        cv::line(topviewImg, cv::Point(0, y), cv::Point(imageSize.width, y), cv::Scalar(255, 0, 0));
    }

    // display image
    string windowName = "3D Objects";
    cv::namedWindow(windowName, 1);
    cv::imshow(windowName, topviewImg);

    if(bWait)
    {
        cv::waitKey(0); // wait for key to be pressed
    }
}


// associate a given bounding box with the keypoints it contains
void clusterKptMatchesWithROI(BoundingBox &boundingBox, std::vector<cv::KeyPoint> &kptsPrev, std::vector<cv::KeyPoint> &kptsCurr, std::vector<cv::DMatch> &kptMatches)
{
    std::vector<double> euclideanDist;
    double distThrehold = 2.0;

    // Find the median euclidean distance between kptsPrev and kptsCurr
    for (auto it = kptMatches.begin(); it != kptMatches.end(); ++it)
    {
        if (boundingBox.roi.contains(kptsCurr[it->trainIdx].pt))
        {
            euclideanDist.push_back(cv::norm(kptsPrev[it->queryIdx].pt - kptsCurr[it->trainIdx].pt));
        }
    }
    const auto med_it = euclideanDist.begin() + euclideanDist.size() / 2;
    const auto med_it_1 = euclideanDist.begin() + euclideanDist.size() / 2 - 1;
    std::nth_element(euclideanDist.begin(), med_it, euclideanDist.end());
    std::nth_element(euclideanDist.begin(), med_it_1, euclideanDist.end());
    double medEuclideanDist = euclideanDist.size() % 2 == 0 ? (*med_it + *med_it_1) / 2.0 : *med_it;

    
    for (auto it = kptMatches.begin(); it != kptMatches.end(); ++it)
    {
        if (boundingBox.roi.contains(kptsCurr[it->trainIdx].pt) &&
            abs(cv::norm(kptsPrev[it->queryIdx].pt - kptsCurr[it->trainIdx].pt) - medEuclideanDist) < distThrehold)
        {
            boundingBox.kptMatches.push_back(*it);
        }
    }
}


// Compute time-to-collision (TTC) based on keypoint correspondences in successive images
void computeTTCCamera(std::vector<cv::KeyPoint> &kptsPrev, std::vector<cv::KeyPoint> &kptsCurr, 
                      std::vector<cv::DMatch> kptMatches, double frameRate, double &TTC, cv::Mat *visImg)
{
    // compute distance ratios between all matched keypoints
    vector<double> distRatios; // stores the distance ratios for all keypoints between curr. and prev. frame
    for (auto it1 = kptMatches.begin(); it1 != kptMatches.end() - 1; ++it1)
    { // outer keypoint loop

        // get current keypoint and its matched partner in the prev. frame
        cv::KeyPoint kpOuterCurr = kptsCurr.at(it1->trainIdx);
        cv::KeyPoint kpOuterPrev = kptsPrev.at(it1->queryIdx);

        for (auto it2 = kptMatches.begin() + 1; it2 != kptMatches.end(); ++it2)
        { // inner keypoint loop

            double minDist = 100.0; // min. required distance

            // get next keypoint and its matched partner in the prev. frame
            cv::KeyPoint kpInnerCurr = kptsCurr.at(it2->trainIdx);
            cv::KeyPoint kpInnerPrev = kptsPrev.at(it2->queryIdx);

            // compute distances and distance ratios
            double distCurr = cv::norm(kpOuterCurr.pt - kpInnerCurr.pt);
            double distPrev = cv::norm(kpOuterPrev.pt - kpInnerPrev.pt);

            if (distPrev > std::numeric_limits<double>::epsilon() && distCurr >= minDist)
            { // avoid division by zero

                double distRatio = distCurr / distPrev;
                distRatios.push_back(distRatio);
            }
        } // eof inner loop over all matched kpts
    }     // eof outer loop over all matched kpts

    // only continue if list of distance ratios is not empty
    if (distRatios.size() == 0)
    {
        TTC = NAN;
        return;
    }

    // compute camera-based TTC from distance ratios
    double medianDistRatio;
    if (distRatios.size() % 2 == 0){
        const auto med_it = distRatios.begin() + distRatios.size() / 2;
        const auto med_it_1 = distRatios.begin() + distRatios.size() / 2 - 1;
        std::nth_element(distRatios.begin(), med_it, distRatios.end());
        std::nth_element(distRatios.begin(), med_it_1, distRatios.end());
        medianDistRatio = (*med_it + *med_it_1) / 2.0;
    } else {
        const auto med_it = distRatios.begin() + distRatios.size() / 2;
        std::nth_element(distRatios.begin(), med_it, distRatios.end());
        medianDistRatio = *med_it;
    }
    


    double dT = 1 / frameRate;
    TTC = -dT / (1 - medianDistRatio);
}


void computeTTCLidar(std::vector<LidarPoint> &lidarPointsPrev,
                     std::vector<LidarPoint> &lidarPointsCurr, double frameRate, double &TTC)
{
    // auxiliary variables
    double dT = 1.0 / frameRate; // time between two measurements in seconds
    double laneWidth = 4.0;      // assumed width of the ego lane

    vector<double> xPrev;

    for (auto it = lidarPointsPrev.begin(); it != lidarPointsPrev.end(); ++it)
    {
        // remove lidar points outside the lane
        if (abs(it->y) <= (laneWidth / 2.0))
        {
            xPrev.push_back(it->x);
        }
    }

    // Since the we want to remove the outlier, we use the 5% percentile of the 
    // x position to represent the closest lidar point to the ego car.
    auto five_per_it = xPrev.begin() + xPrev.size() / 20;
    std::nth_element(xPrev.begin(), five_per_it, xPrev.end());
    double five_per_XPrev = *five_per_it;


    vector<double> xCur;

    for (auto it = lidarPointsCurr.begin(); it != lidarPointsCurr.end(); ++it)
    {
        // remove lidar points outside the lane
        if (abs(it->y) <= (laneWidth / 2.0))
        {
            xCur.push_back(it->x);
        }
    }

    // Since the we want to remove the outlier, we use the 5% percentile of the 
    // x position to represent the closest lidar point to the ego car.
    five_per_it = xCur.begin() + xCur.size() / 20;
    std::nth_element(xCur.begin(), five_per_it, xCur.end());
    double five_per_XCur = *five_per_it;

    // compute TTC from both measurements
    if (abs(five_per_XPrev - five_per_XCur) > std::numeric_limits<double>::epsilon())
    {
        TTC = five_per_XCur * dT / (five_per_XPrev - five_per_XCur);
    }
    else
    {
        TTC = NAN;
    }    
}


void matchBoundingBoxes(std::vector<cv::DMatch> &matches, std::map<int, int> &bbBestMatches, DataFrame &prevFrame, DataFrame &currFrame)
{
    vector<int> matchedBoxIDs;
    
    // Loop over all bounding boxed in the current frame.
    for (auto cur_it = currFrame.boundingBoxes.begin(); cur_it != currFrame.boundingBoxes.end(); ++cur_it)
    {
        // Stores match count, boxID -> count.
        map<int, int> bbMatchCount;

        // Loop over all matches.
        for (auto match_it = matches.begin(); match_it != matches.end(); ++match_it)
        {
            // Current frame contains match point
            if (cur_it->roi.contains(currFrame.keypoints[match_it->trainIdx].pt))
            {
                // Loop over bounding boxes in the previous frame
                for (auto prev_it = prevFrame.boundingBoxes.begin(); prev_it != prevFrame.boundingBoxes.end(); ++ prev_it)
                {
                    // Update the bbMatchCount if the previous bounding box also contains the keypoint.
                    if (prev_it->roi.contains(prevFrame.keypoints[match_it->queryIdx].pt))
                    {
                        if (bbMatchCount.count(prev_it->boxID) == 0)
                        {
                            bbMatchCount.insert(pair<int, int>(prev_it->boxID, 0));
                        }
                        bbMatchCount.at(prev_it->boxID)++;
                    }
                }
            }
        }

        int bestPrevBoxID;
        int maxMatchCount = 0;
        for (auto count_it = bbMatchCount.begin(); count_it != bbMatchCount.end(); ++count_it)
        {
            if ((count_it->second > maxMatchCount) && (find(matchedBoxIDs.begin(), matchedBoxIDs.end(), count_it->first) == matchedBoxIDs.end()))
            {
                maxMatchCount = count_it->second;
                bestPrevBoxID = count_it->first;
            }
        }
        if (maxMatchCount > 0)
        {
            matchedBoxIDs.push_back(bestPrevBoxID);
            bbBestMatches.insert(pair<int, int>(bestPrevBoxID, cur_it->boxID));
        }


    }
}