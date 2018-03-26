This package is a filter that selects events from the Mu primary dataset from the CMS open
data release. In particular, the event is selected if (for simplicity) there are precisely two muons
in the event and at least one is a global muon. A csv file containing the four-vector information, charge, and invariant mass of the two muons is produced.

See http://opendata.cern.ch for more information and for context on the instructions below.

To produce files in the VM open a terminal with the X terminal emulator (an icon bottom-left of the VM screen)
and input the commands as explained below.

* Create a CMSSW environment: 

```
    cmsrel CMSSW_5_3_32
```

* Change to the CMSSW_5_3_32/src/ directory:

```
    cd CMSSW_5_3_32/src/
```
* Initialize the CMSSW environment:

```
    cmsenv
```
* Clone the source code:

```
    git clone git://github.com/tpmccauley/dimuon-filter DimuonFilter/DimuonFilter
```

* Checkout the 2011A branch:

```
    cd DimuonFilter/DimuonFilter
    git checkout 2011A
```

* Compile the code with the command:

```
    scram b
```

* Run the example configuration file (see comments in the file on changing parameters):

```
    cmsRun DoubleMuRun2011A.py
```
which will produce a csv file.
Enjoy!
