# Paper
> F. Liu, C. Masouros, A. Li, H. Sun and L. Hanzo, "MU-MIMO Communications With MIMO Radar: From Co-Existence to Joint Transmission," in IEEE Transactions on Wireless Communications, vol. 17, no. 4, pp. 2755-2770, April 2018. 

- All codes are contributed by Nate Raymondi (Graduate Student, Rice University), 8/13/2020
- [Click here to preview this paper in IEEE](https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=8288677)


# Abstract
Beamforming techniques are proposed for a joint multi-input-multi-output (MIMO) radar-communication (RadCom) system, where a single device acts as radar and a communication base station (BS) by simultaneously communicating with downlink users and detecting radar targets. Two operational options are considered, where we first split the antennas into two groups, one for radar and the other for communication. Under this deployment, the radar signal is designed to fall into the nullspace of the downlink channel. The communication beamformer is optimized such that the beampattern obtained matches the radar s beampattern while satisfying the communication performance requirements. To reduce the optimizations  constraints, we consider a second operational option, where all the antennas transmit a joint waveform that is shared by both radar and communications. In this case, we formulate an appropriate probing beampattern, while guaranteeing the performance of the downlink communications. By incorporating the SINR constraints into objective functions as penalty terms, we further simplify the original beamforming designs to weighted optimizations, and solve them by efficient manifold algorithms. Numerical results show that the shared deployment outperforms the separated case significantly, and the proposed weighted optimizations achieve a similar performance to the original optimizations, despite their significantly lower computational complexity.

# About this code
### Software platform
- MATLAB 2014/2016/2018/2020
- To run those codes, please download and install [CVX](http://cvxr.com/cvx/)

### Content
- To product Fig.3-Fig.4 in this paper

### Usage
1. Due to the version of your MATLAB, variable name "sigma" is the name of an existing MATLAB function or directory. You need choose a different name (sigma_temp for example)
2. Run `SeparatedDeployment_BPmatching.m`
3. Run `SeparatedDeployment_SidelobeMin.m`
