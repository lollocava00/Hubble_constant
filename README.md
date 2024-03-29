# Hubble_constant
Determining the Value for the Hubble constant from a sample of Supernovae light curves given by the Canergie supernova project 2011.
Starting from the available data(70 SNe .dat files), only for supernovae having data before the maximum, calculate:

• the time of the maximum,

• the apparent magnitude in b-band at the maximum,

• the b-band magnitude 15 days after the maximum.

From these data we then obtain:

• the absolute magnitude in band B, correcting for extinction,

• the module of distance μ,

• the distance (brightness) of the object.

Finally, use the distance-redshift relationship to estimate the constant of Hubble H0.

The full description of the code and the numerical methods used is in italian in 'description.pdf' file.

Some examples of Sne light curves:

![Alt text](Results/grafico1.png?raw=true)

![Alt text](Results/grafico5.png?raw=true)

Example of Montecarlo Technique on 3 Sne light curves:

![Alt text](Results/grafico3.png?raw=true)

The distance modulus-redshift relation:

![Alt text](Results/grafico2.png?raw=true)
