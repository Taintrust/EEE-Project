function mpc = case15nbr
%CASE15NBR  Power flow data for 15 bus distribution system from Battu, et al
%   Please see CASEFORMAT for details on the case file format.
%
%   Data from ...
%       Battu NR, Abhyankar AR, Senroy N (2016) DG Planning with Amalgamation
%       of Operational and Reliability Considerations. Int J Emerg Electr
%       Power Syst 17:131-141. doi: 10.1515/ijeeps-2015-0142
%       URL: https://doi.org/10.1515/ijeeps-2015-0142

%% MATPOWER Case Format : Version 2
mpc.version = '2';

%%-----  Power Flow Data  -----%%
%% system MVA base
mpc.baseMVA = 100;

%% bus data
%%type 1  = PQ bus, type 2 = PV bus and type 3 = slack (reference bus)
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
mpc.bus = [ %% (Pd and Qd are specified in kW & kVAr here, converted to MW & MVAr below)
	1	3	0	0	0	0	1	1	0	11	1	1.2	0.8; % Vmax and Vmin updated
	2   2   44.1    44.991  0   0   1   1   0   11  1   1.2 0.8; % Changed to PV bus
	3	1	70	71.4143	0	0	1	1	0	11	1	1.2	0.8;
	4	1	5380	142.8286	0	0	1	1	0	11	1	1.2	0.8; % Node 4 changed to PV bus 5.2MW + 140
	5	1	44.1	44.991	0	0	1	1	0	11	1	1.2	0.8;
	6	1	140	142.8286	0	0	1	1	0	11	1	1.2	0.8;
	7	1	140	142.8286	0	0	1	1	0	11	1	1.2	0.8;
	8	1	70	71.4143	0	0	1	1	0	11	1	1.2	0.8;
	9	1	70	71.4143	0	0	1	1	0	11	1	1.2	0.8;
	10	1	44.1	44.991	0	0	1	1	0	11	1	1.2	0.8;
	11	1	140	142.8286	0	0	1	1	0	11	1	1.2	0.8;
	12	1	70	71.4143	0	0	1	1	0	11	1	1.2	0.8;
	13	1	44.1	44.991	0	0	1	1	0	11	1	1.2	0.8;
	14	1	70	71.4143	0	0	1	1	0	11	1	1.2	0.8;
	15	1	140	142.8286	0	0	1	1	0	11	1	1.2	0.8;
];

%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	Pc1	Pc2	Qc1min	Qc1max	Qc2min	Qc2max	ramp_agc	ramp_10	ramp_30	ramp_q	apf
mpc.gen = [
	1	0	0	10	-10	1	100	1	10	0	0	0	0	0	0	0	0	0	0	0	0;  
	2   1. 49   50  -50  1   100   1   10     0   0   0   0   0   0   0   0   0   0   0   0; % New CHP generator at node 2, 1.488MW
];

%% branch data
%	fbus	tbus	r	x	b	rateA	rateB	rateC	ratio	angle	status	angmin	angmax
mpc.branch = [
    1	2	0.7766	0.7596	0	700	0	0	0	0	1	-360	360;    %%Updated rateA to 700MW
	2	3	0.6716	0.6569	0	700	0	0	0	0	1	-360	360;    %%rateA is maximum power flow capacity of a branch
	3	4	0.4827	0.4722	0	700	0	0	0	0	1	-360	360;
	4	5	0.8744	0.5898	0	700	0	0	0	0	1	-360	360;
	2	9	1.1554	0.7793	0	700	0	0	0	0	1	-360	360;
	9	10	0.9680	0.6529	0	700	0	0	0	0	1	-360	360;
	2	6	1.4677	0.9900	0	700	0	0	0	0	1	-360	360;
	6	7	0.6245	0.4213	0	700	0	0	0	0	1	-360	360;
	6	8	0.7182	0.4844	0	700	0	0	0	0	1	-360	360;
	3	11	1.0305	0.6951	0	700	0	0	0	0	1	-360	360;
	11	12	1.4052	0.9478	0	700	0	0	0	0	1	-360	360;
	12	13	1.1554	0.7793	0	700	0	0	0	0	1	-360	360;
	4	14	1.2803	0.8636	0	700	0	0	0	0	1	-360	360;
	4	15	0.6870	0.4634	0	700	0	0	0	0	1	-360	360;
];



%%-----  OPF Data  -----%%
%% generator cost data
%	1	startup	shutdown	n	x1	y1	...	xn	yn
%	2	startup	shutdown	n	c(n-1)	...	c0
mpc.gencost = [
	2	0	0	3	0	20	0;
	2   0   0   3   0.01   15   0; % Cost for CHP generator at Node 2
];


%% convert branch impedances from Ohms to p.u.
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;

%% convert loads from kW to MW
mpc.bus(:, [PD, QD]) = mpc.bus(:, [PD, QD]) / 1e3;