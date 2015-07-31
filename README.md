# Access to longitudinal HIVEVO data

This is a simple collection of classes that provide a python interface to precomputed data from the HIVEVO project.

## Configuration
The repository does not contain any data itself and relies on an external data source.

The location of the root data folder can be specified via an environment variable, `HIVEVO_ROOT_DATA_FOLDER`. It can be set either in the shell (e.g. in .bashrc/.bash_profile) or after python startup, via:
```python
import os
os.environ['HIVEVO_ROOT_DATA_FOLDER'] = '/your/location/on/disk/'
```

## Patient
To load a patient:

```Python
from hivevo.patients import Patient; patient = Patient.load('p1')
```

## Sample



