# API reference

## Phases

```@docs
cellPhases
PhaseResult
Track
bandpassFilter
timeShift
```

## Readouts

```@docs
cycleLengths
CycleLengths
phaseValueCheck
PhaseValueCheck
phaseCoherence
PhaseCoherence
orderParameter
OrderParameter
phaseWaves
PhaseWaves
cleanShare
```

## Replicates

```@docs
Replicates
replicateSpread
```

## Plot data

```@docs
roseData
coherenceData
orderParameterData
cycleLengthData
phaseValueData
traceData
portraitData
```

## Makie recipes

Loaded with a Makie backend (`using CairoMakie` or similar).

```@docs
roseplot
roseplot!
coherenceplot
coherenceplot!
orderparameterplot
orderparameterplot!
cyclelengthplot
cyclelengthplot!
phasevalueplot
phasevalueplot!
portraitplot
portraitplot!
```

## Synthetic data

```@docs
agentField
planeWave
pulseTrain
oscillatorField
```

## PhysiCell

```@docs
cellPhases(::PhysiCellOutput.PhysiCellSequence)
```

## Module

```@docs
CellPhases
```
