

set Generators; 		#indexed by g
set Buses;				#indexed by b
set Lines;				#indexed by l

param Cost1{g in Generators};
param Cost2{g in Generators};
param Cost3{g in Generators};
param RealPowerGenLimitLower{g in Generators};
param RealPowerGenLimitUpper{g in Generators};
param ReactivePowerGenLimitLower{g in Generators};
param ReactivePowerGenLimitUpper{g in Generators};
param GenLocation{g in Generators}; 

param SusceptanceOfReactivePower{b in Buses};
param VoltageLimitUpper{b in Buses};
param VoltageLimitLower{b in Buses};
param GeneratorIndex{b in Buses};
																				
param RealPowerBusDemand{b in Buses};
param ReactivePowerBusDemand{b in Buses};

param ThermalLineLimit{l in Lines};

param LineConductance{l in Lines};
param LineSusceptance{l in Lines};
param LineResistance{l in Lines};
param ShuntSusceptance{l in Lines};

param FromBus{l in Lines};
param ToBus{l in Lines};

# Upper an lower limits defined here in the declaration for
# the Voltage Level, for each bus, and for Real and Reactive
# Power Generation, for each generator
var VoltageLevel{b in Buses} >= VoltageLimitLower[b], <= VoltageLimitUpper[b];
var VoltagePhase{b in Buses};
var RealPowerGeneration{g in Generators} <= RealPowerGenLimitUpper[g],>= RealPowerGenLimitLower[g];
var ReactivePowerGeneration{g in Generators} <= ReactivePowerGenLimitUpper[g],>= ReactivePowerGenLimitLower[g];
var RealPowerGenerationAtBuses{b in Buses};
var ReactivePowerGenerationAtBuses{b in Buses};
var RealPowerInjectionFromBus{l in Lines};		
var RealPowerInjectionToBus{l in Lines};				
var ReactivePowerInjectionFromBus{l in Lines};	
var ReactivePowerInjectionToBus{l in Lines};				
														

#############################

# Objective function
minimize OperatingCost: sum{g in Generators}( Cost1[g]*RealPowerGeneration[g]^2 + Cost2[g]*RealPowerGeneration[g] + Cost3[g]);

#############################


# Kirchhoff Current Law constraints
# 'AtBuses' constraints are just sums  to make the KCL constraints more readable
subject to RealPowerGenerationAtBusesCS{b in Buses}:
	RealPowerGenerationAtBuses[b] = sum{g in Generators: GenLocation[g]==b} RealPowerGeneration[g];

subject to ReactivePowerGenerationAtBusesCS{b in Buses}:
	ReactivePowerGenerationAtBuses[b] = sum{g in Generators: GenLocation[g]==b} ReactivePowerGeneration[g];

subject to RealCurrentCS{b in Buses}:
	RealPowerGenerationAtBuses[b] = RealPowerBusDemand[b] + sum{l in Lines: FromBus[l]==b} RealPowerInjectionFromBus[l] + sum{l in Lines: ToBus[l] == b} RealPowerInjectionToBus[l];

subject to ReactiveCurrentCS{b in Buses}:
	ReactivePowerGenerationAtBuses[b] = ReactivePowerBusDemand[b]  + sum{l in Lines: FromBus[l]==b} ReactivePowerInjectionFromBus[l] + sum{l in Lines: ToBus[l] == b} ReactivePowerInjectionToBus[l] + SusceptanceOfReactivePower[b]*VoltageLevel[b]^2; 
	
# Kirchoff Voltage Law constraints
# Defined at both ends of the line
subject to RealVoltageCS1{l in Lines}:
   RealPowerInjectionFromBus[l] = LineConductance[l]*(VoltageLevel[FromBus[l]])^2-VoltageLevel[FromBus[l]]*VoltageLevel[ToBus[l]]*(LineConductance[l]*cos(VoltagePhase[FromBus[l]]-VoltagePhase[ToBus[l]])+LineSusceptance[l]*sin(VoltagePhase[FromBus[l]]-VoltagePhase[ToBus[l]]));

subject to RealVoltageCS2{l in Lines}:
   RealPowerInjectionToBus[l] = LineConductance[l]*(VoltageLevel[ToBus[l]])^2-VoltageLevel[ToBus[l]]*VoltageLevel[FromBus[l]]*(LineConductance[l]*cos(VoltagePhase[ToBus[l]]-VoltagePhase[FromBus[l]])+LineSusceptance[l]*sin(VoltagePhase[ToBus[l]]-VoltagePhase[FromBus[l]]));   

subject to ReactiveVoltageCS1{l in Lines}: 
	ReactivePowerInjectionFromBus[l] = - (LineSusceptance[l]+ShuntSusceptance[l]/2)*(VoltageLevel[FromBus[l]])^2-VoltageLevel[FromBus[l]]*VoltageLevel[ToBus[l]]*(LineConductance[l]*sin(VoltagePhase[FromBus[l]]-VoltagePhase[ToBus[l]])-LineSusceptance[l]*cos(VoltagePhase[FromBus[l]]-VoltagePhase[ToBus[l]]));

subject to ReactiveVoltageCS2{l in Lines}: 
	ReactivePowerInjectionToBus[l] = - (LineSusceptance[l]+ShuntSusceptance[l]/2)*(VoltageLevel[ToBus[l]])^2-VoltageLevel[ToBus[l]]*VoltageLevel[FromBus[l]]*(LineConductance[l]*sin(VoltagePhase[ToBus[l]]-VoltagePhase[FromBus[l]])-LineSusceptance[l]*cos(VoltagePhase[ToBus[l]]-VoltagePhase[FromBus[l]]));
		
# Thermal line limits constraints
# Defined at both ends of the lines
subject to ThermalLimitsCS1{l in Lines}:
	(ThermalLineLimit[l])^2 >= (RealPowerInjectionFromBus[l])^2 + (ReactivePowerInjectionFromBus[l])^2;

subject to ThermalLimitsCS2{l in Lines}:
	(ThermalLineLimit[l])^2 >= (RealPowerInjectionToBus[l])^2 + (ReactivePowerInjectionToBus[l])^2;
	
# Set voltage phase at a slack bus equal to 0 
# We chose to set the location of generator 1 as the slack bus
# (Arbitrary choice of a generation bus)
#subject to VoltagephaseCS{g in Generators: g == 1}:
#	VoltagePhase[GenLocation[g]] == 0;
subject to VoltagePhaseCS{b in Buses: b == 1}:
	VoltagePhase[b] == 0;	