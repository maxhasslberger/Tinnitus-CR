# Tinnitus_CR

Generate individualized Acoustic Coordinated Reset (aCR) wav files based on [Tass et al., 2019](https://www.nature.com/articles/s41598-019-49945-w) and [Munjal et al., 2021](https://www.frontiersin.org/articles/10.3389/fnetp.2021.734344/full).

`tinnitus_cr.m`

1. Adapt `f_t` based on your Tinnitus center frequency
2. (Adapt other init parameters arbitrarily)
3. Run script and obtain individual wav file

## Modes
Both regular and noisy aCR protocols are supported. Additionally, one can pick the `log` mode to logarithmically assign stimulus frequencies around the Tinnitus frequency.
