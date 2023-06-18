# MopacWIthLibra
A place to keep and update code


	\section{Non-adiabatic Molecular Dynamics for INDO}
	For this study we implemented the necessary post-processing methods to calculate non-adiabatic couplings in the molecular orbital, Slater determinant(SD), and configuration interaction(CI) bases.
	We focus on doing non-adiabatic molecular dynamics(NA-MD) where the nuclei are treated classically and the adiabatic electronic states are parametrically dependent on the atomic coordinates which vary with time.
	The adibatic electronic states or molecular orbitals are defined as a linear combination of atomic orbitals
	\begin{equation}
		\phi_i(r,R,t)=\sum_{i} c_i \chi_i(r;R(t)) 
	\end{equation}
	which are dependent on the electron position (r) and parametrically dependent on the classical trajectory(R(t)) of the nuclei.
	The non-adiabatic coupling which measures the interaction between electronic and nuclear degrees of freedom.
	We calculate the time-derivative non-adiabatic coupling ($d_{ij}$)(TNAC) as defined by Hammes-Schiffer and Tully\autocite{hammesschiffer_proton_1994}
	\begin{equation}
		d_{ij}(t+\frac{dt}{2}) = -i \hbar \frac{\langle \phi_i(t) | \phi_j(t+dt)\rangle - \langle \phi_i(t+dt) | \phi_j(t) \rangle }{2\cdot dt } \label{eq:full_NAC}.
	\end{equation}
	The non-adiabatic couplings are then used to build the vibrational Hamiltonian, $H_{vib}$ which contains the non-adiabatic couplings and average energy of the states between time steps.
	\begin{equation}
		H_{vib(ii)} = \frac{E_i(t) + E_i(t+dt)}{2}
		\label{eq:Hvib_off}
	\end{equation}
	\begin{equation}
		H_{vib(ij)} = d_{ij}
		\label{eq:Hvib_dia}
	\end{equation}
	
	
	The process to build the vibrational Hamiltonian matrix in the molecular orbital basis is shown above.
	However, when calculating TNACs in the SD basis one must consider the overlap between all of the occupied orbitals between two SDs.
	A SD is defined as
	\begin{equation}
		\Phi_p = \frac{1}{\sqrt{N!}} det
		\begin{vmatrix}
			\phi_{p_1}(1,t) & \dots & \phi_{p_1}(N,t) \\
			\vdots & \ddots & \vdots \\
			\phi_{p_i}(1,t) & \dots & \phi_{p_i}(N,t)
		\end{vmatrix}
		\label{eq:SD}.
	\end{equation}
	The SD is a convenient way to represent a collection of one electron wavefunctions.
	The indexing in $\Phi_p$ runs from 1 to i, where i is the number of molecular orbital coefficients.
	There are N electrons represented by this SD and these one electron wavefunctions are dependent on the time steps in the classical trajectory.
	
	
	To evaluate the overlap between SDs\autocite{plasser_efficient_2016} $\langle \Phi_p|\Phi_q \rangle$ a matrix that consist of overlaps between all of the molecular orbitals in $\Phi_p$ and $\Phi_q$ is constructed and the determinant of that matrix yields the overlap between SDs \ref{eq:SD_ovlp}. 
	\begin{equation}
		\langle \Phi_p | \Phi_q \rangle =  \begin{vmatrix}
			\langle \phi_{p_1}(t)| \phi_{q_1}(t)\rangle & \dots & \langle \phi_{p_1}(t)| \phi_{q_n}(t) \rangle \\
			\vdots & \ddots & \vdots \\
			\langle \phi_{p_n}(t)|\phi_{q_1}(t) \rangle & \dots & \langle \phi_{p_n}(t) | \phi_{q_n}(t) \rangle \end{vmatrix}
		\label{eq:SD_ovlp}
	\end{equation}
	Here there are n molecular orbitals which also depend on the classical trajectory.
	The time overlaps($S_t$) in the SD basis can be computed in a similar manner by substituting $\Phi_q$ with $\Phi_q^\prime$($\Phi_q$(t+dt)) and using the same method in eq. \ref{eq:SD_ovlp}.
	
	
	Expanding to the configuration interaction basis involves a linear transformation of the time overlap matrices in the SD basis to the CI basis(eq.\ref{eq:CI_ovlp}).
	\begin{equation}
		S_{ij}(t,t+dt) =\langle \Psi_i(t)|\Psi_j(t+\Delta t)\rangle = T^\intercal(t)\langle \Phi_{I=i_1,i_2,\cdots,i_n}(t) | \Phi_{J=j_1,j_2,\cdots,j_n}(t+\Delta t) \rangle T(t+\Delta t) \label{eq:CI_ovlp}
	\end{equation}
	Here $\langle \Psi_i(t)|\Psi_j(t+\Delta t)\rangle$ is the time overlap between CI states i and j and each of these CI states involves a set of SDs($J=j_1,j_2,\cdots,j_n$).
	The number of SDs computed limits the number of excited states that can be computed therefore the size of the T vector is limited by the number of excited states.
	The T vector contains the CI coefficients for all SDs that are used to describe an excited state.
	From the time overlaps we can calculate the TNACs. 
	\begin{equation}
		d_{ij} = \frac{S_{ij}(t,t+dt) - S_{ji}^\ast(t,t+dt)}{2\Delta t}.
		\label{eq:TNACs_S}
	\end{equation}
	
	The vibrational Hamiltonian will have the same form as in equations \ref{eq:Hvib_off} and \ref{eq:Hvib_dia}, but depending on the basis, the energies and operations to compute the overlaps will be different.
	All code necessary to extract energies and wavefunctions from Mopac output files and calculate the overlaps and TNACs has been written.
	For this study a comparison between the results produced by INDO and DFT still needs to done for this project.
	Moving forward with this project we will be able to study photophysical processes of molecules and metal clusters using INDO and the FSSH algorithm to study the decay of plasmon-like excited.
	In addition to studying the non-radiative decay of plasmonic-like states our INDO calculations will take a fraction of the time compared to DFT.


