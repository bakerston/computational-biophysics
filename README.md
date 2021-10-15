## hbond_data




##### Read pdb file and calculate time dependence hydrogen bonds between any two amino chains.

`./hbond_data.py file_name.pdb`



#### Read pdb file and calculate time dependence contact number within C-Alpha between any two amino chains.
`./hbond_ca_con.py`



#### Read Secondary-Structure output file and calculate averaged SS.

`./dssp.py`
|list|type|structure|Amnio Acids|
|---|---|---|---|
|`alist`|a alpha-helix|H, G, I|His, Gly, Ile|
|`blist`|b beta-sheet|B, E|Glu, Asx|
|`clist`|t tertiary|T|Thr|
|`dlist`|c coil|S, C|Ser, Cys| 


![](https://github.com/bakerston/computational-biophysics/blob/master/src/hbonds.png)
![](https://github.com/bakerston/computational-biophysics/blob/master/src/hbond_index.png)
#### Read pdb file and calculate time dependence Q value of the selected chain.
`q_value_hbond.py `

| inputs | parameters | default value |
|---|---|---|
dis_cut| cut off distance |3.5 angstrom
ang_cut| cut off angle |120 degree
res_num| residue number on the chain|37 N/A
eff_range| effective range |0 angstrom

![](https://github.com/bakerston/computational-biophysics/blob/master/src/ss.png)

#### Create data for (Origin Software)
`xyz_grid.py`


#### Elongate protein fibrils by adding chains and maintain double helix structure.
`firbil_elong.py`
![](https://github.com/bakerston/computational-biophysics/blob/master/src/elongation.png)
![]()
