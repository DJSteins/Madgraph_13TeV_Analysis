#!/bin/bash

currentDir=$PWD

cd /home/dsteiner/madgraph/MG5_aMC_v2_6_3_2/Delphes

analysis_path=/home/dsteiner/madgraph/analysis

sr=tN_diag_high
run=01
csvFile="$sr"_"TTBAR_1Lepton.csv"

root -b<<EOF
gSystem->Load("libDelphes");
gInterpreter->AddIncludePath("/home/dsteiner/RestFrames/inc")
gInterpreter->AddIncludePath("/home/dsteiner/gambit/Utils/include")
gInterpreter->AddIncludePath("/home/dsteiner/gambit/Logs/include")
gInterpreter->AddIncludePath("/home/dsteiner/gambit/contrib/heputils/include")
gInterpreter->AddIncludePath("/home/dsteiner/gambit/contrib/mcutils/include")
gInterpreter->AddIncludePath("/home/dsteiner/gambit/ColliderBit/include/")
gInterpreter->AddIncludePath("/usr/lib/gcc/x86_64-redhat-linux/4.4.7/include/")
gInterpreter->LoadFile("/home/dsteiner/gambit/ColliderBit/src/analyses/mt2_bisect.cpp")
gInterpreter->LoadFile("/home/dsteiner/gambit/ColliderBit/src/Utils.cpp")
gInterpreter->Load("/usr/lib/gcc/x86_64-redhat-linux/4.4.7/libgomp.so")
.X $analysis_path/tev13_analysis.C("/home/dsteiner/madgraph/MG5_aMC_v2_6_3_2/$sr/Events/run_$run/tag_1_delphes_events.root","$analysis_path/$csvFile")
.q
EOF

cd $currentDir
