# Paper
> Fan Liu*; Longfei Zhou; Christos Masouros; Ang Li; Wu Luo; Athina Petropulu: Toward Dual-functional Radar-Communication Systems: Optimal Waveform Design, IEEE Transactions on Signal Processing, vol. 66, no. 16, pp. 4264-4279, Aug. 2018. (ranked 8th most popular article in IEEE TSP, Feb 2021)

- All codes are contributed by Fan Liu, SUSTech
- [Click here to preview this paper in IEEE](https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=8386661)


# Abstract
We focus on a dual-functional multi-input-multioutput (MIMO) radar-communication (RadCom) system, where a single transmitter with multiple antennas communicates with downlink cellular users and detects radar targets simultaneously. Several design criteria are considered forminimizing the downlink multiuser interference. First, we consider both omnidirectional and directional beampattern design problems, where the closedform globally optimal solutions are obtained. Based on the derived waveforms, we further consider weighted optimizations targeting a flexible tradeoff between radar and communications performance and introduce low-complexity algorithms. Moreover, to address the more practical constant modulus waveform design problem, we propose a branch-and-bound algorithm that obtains a globally optimal solution, and derive its worst-case complexity as function of the maximum iteration number. Finally, we assess the effectiveness of the proposed waveform design approaches via numerical results.

# About this code
### Software platform
- MATLAB 2016/2018/2020
- To run those codes, please download and install [CVX](http://cvxr.com/cvx/) & [Manopt](https://www.manopt.org/)
### Content
- The folder `Waveform Design With Given Radar Beampatterns` is provided to produuct Fig.3-Fig.7
- The folder `Constant modulus` is provided to product Fig.8-Fig.10
