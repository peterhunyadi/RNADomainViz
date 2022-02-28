# RNADomainViz: Drawing domain-level nucleic acid structures
Software Development Course Project 2021W

Domain-level nucleic acid structures can be visualized by many tools already or by hand. This tool can render 2D representations of  nucleic acid molecules from kernel string input with extra layers of precision.
This extra inputs are :
forcing angles between domains with in-line additions
specifying domain length from outside of the string

You can also choose a color scheme from 4 pre-added colorblind friendly palettes (“IBM”, “Wong”, “Magma”, “Plasma”.

##Packages
The following packages are called: numpy, networkx, shapely, drawSvg

## How to use:
Just import the domain_visualization function from domainviz

```python

from domainviz import domain_visualization

string = " aa( i1 da( dl ) i2 aca( acl ) vr pa( pl ) i3 ) c"
filename = "trna1.svg"

master_function(string, length_input_dict, palette, filename, orient=-1)
```

Can force some angles and length like that:
```python

string = "@-90 aa( i1 da( dl ) i2 aca( acl ) vr pa( pl ) i3 ) @0 c"
length_input_dict = {0:30, 2:20, 6:20, 10:25}
filename = "trna2.svg"

master_function(string, length_input_dict=length_input_dict, filename=filename)
```

Can change the other options, orientation or palette:
```python

string = "@-90 aa( i1 da( dl ) i2 aca( acl ) vr pa( pl ) i3 ) @0 c"
length_input_dict = {0:30, 2:20, 6:20, 10:25}
palette = "Plasma"
filename = "trna.svg"

master_function(string, length_input_dict=length_input_dict, palette = palette, filename=filename, orient=-1)
```
