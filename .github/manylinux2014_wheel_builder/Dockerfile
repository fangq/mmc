FROM quay.io/pypa/manylinux2014_x86_64
RUN yum install -y ocl-icd ocl-icd-devel vim-common && yum clean all
ENV PATH /opt/rh/devtoolset-7/root/usr/bin:${PATH}
COPY build_wheels.sh /
WORKDIR /src
ENTRYPOINT ["/build_wheels.sh"]
