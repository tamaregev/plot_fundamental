# Matlab library of functions and scripts used for the prosody-meaning project
A behavioral online study testing the detection of pitch changes in spoken phrases, in locations that are relevant (versus not relevant) to meaning.

- Download the folders myspectrogram and STRAIGHT and add them to the matlab path
- Fix the paths inside the various functions and scripts to point to the relevant locations on your machine
    
## Specific functions:
  - make_stim_prosody_meaning.m:
  Pipeline of stimulus creation

  - plot_fundamental.m : 
  Analyze and plot fundamental frequency (f0) of a speech sentence. Instructions:
    - input parameters: start from the values suggested in the function commented initial lines, and change if necessary.
    
  - plotf0morphs.m :
  Analyze and plot f0 of several sounds on top of each other
    
  - logspace4morph.m :
  Finds logarithmically equally spaced frequencies in a given range. Find the morph numbers (out of a linear frequency space, as in STRAIGHT morphs) to use to approximate these. 
  
  - SynthPureTone.m
  - SynthHarmonicTone.m
  Rough functions to synthesize pure and harmonic tones


