# Tasrkleat

Note: `linuxbrew` is not used in Dockerfile because of the
[uid-mapping problem](https://github.com/docker/docker/issues/7198) between
inside and outside a Docker container. Besides, The author finds it inconvenient
to configure specific versions for different bioinformatics softwares with
`linuxbrew`.
