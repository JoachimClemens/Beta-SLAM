# β-SLAM

Simultaneous localization and grid mapping with beta distributions

### Authors

Joachim Clemens, Tobias Kluth, Thomas Reineking


### Description

Simultaneous localization and mapping (SLAM) is one of the most frequently studied problems in mobile robotics.
Different map representations have been proposed in the past and a popular one are occupancy grid maps, which are particularly well suited for navigation tasks. 
The uncertainty in these maps is usually modeled as a single Bernoulli distribution per grid cell.
This has the disadvantage that one cannot distinguish between uncertainty caused by different phenomena like missing or conflicting information. 
In the β-SLAM algorithm, we overcome this limitation by modeling the occupancy probabilities as random variables.
Those are assumed to be beta-distributed and account for the different causes of uncertainty. 
The additional information provided by this approach can be utilized for navigation tasks like path planning or active exploration.


### Paper Describing the Approach

Joachim Clemens, Tobias Kluth, Thomas Reineking, *β-SLAM: Simultaneous localization and grid mapping with beta distributions*, Information Fusion, volume 52, 2019, pages 62-75, [doi:10.1016/j.inffus.2018.11.005](https://doi.org/10.1016/j.inffus.2018.11.005).

### Example Maps

Resulting maps for the Intel Research Lab (raw dataset recorded by Dirk Hähnel) and Cartesium building, University of Bremen (raw dataset recorded by Cyrill Stachniss):

<image src="/images/intel.png" alt="Intel" height="200px" />
<image src="/images/cartesium.png" alt="Cartesium" height="200px" />

The color coding is red for occupied, green for free, blue for the superset (corresponding to unknown areas), and black for empty set (corresponding to conflicts).
The estimated path is shown in yellow.

### Software Requirements

The software is developed and tested on Ubuntu Linux 16.04.
It should run on other recent Linux and UNIX systems as well, while some modifications may be required to run it on Windows.
Furthermore, the following software is required:

* CMake (package `cmake`)
* Qt5 (package `qtbase5-dev`)
* Eigen3 (package `libeigen3-dev`)

### Compilation and Running

Once all required software is installed, the code can be compiled and executed as follows:

```
git clone https://github.com/JoachimClemens/Beta-SLAM.git
cd Beta-SLAM
mkdir build
cd build
cmake ..
make
cd ../bin
./BSlamCarmenGui --filename <carmen-log-file>
```

The input file has to be in CARMEN log format.
Be sure that you configure the sensor parameters with `--laser-start-angle` (default: -90.0) and `--laser-angular-res` (default: 1.0) correctly, which are both given in degrees.
Datasets are, e.g., provided by [Cyrill Stachniss](http://www2.informatik.uni-freiburg.de/~stachnis/datasets.html), while one of the most popular datasets might be the one of the Intel Research Lab recorded by Dirk Hähnel (http://www2.informatik.uni-freiburg.de/~stachnis/datasets/datasets/intel-lab/intel.log.gz, see above).


### License Information

β-SLAM is published under the BSD License. See [LICENSE](LICENSE) for further information.
