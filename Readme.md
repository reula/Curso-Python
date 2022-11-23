## Material de Curso Python en FaMAF

Aquí podrá encontrar el material del curso.

Se está usando material un poco obsoleto pues FeniCs ha cambiado a FeniCsX.
Para hacer funcionar esta versión hay que instalar el FeniCs legacy usando versiones viejas de todo. 
Aún así hay algunas cosas relacionadas a matplotlib que no funcionan correctamente. 

Una vez instalado conda:

```
> conda create -n fenicsproject_old -c conda-forge python=3.8 jupyter fenics
> conda activate fenicsproject_old 
> conda install matplotlib  
```

Hay que instalar por separado el generador de grillas:

```
> conda install -c conda-forge mshr
```
