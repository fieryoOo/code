C######################################################################
        function PHA (y,x)
c----------------------------------------------------------------------
        pi=3.14159
        a=sqrt(x**2+y**2)
        if(a.lt.1.e-08) a=1.e-08
        pha=asin(y/a)
        if(x.lt.0.)   pha=pi-pha
        if(pha.lt.0.) pha=pha+2.*pi
        RETURN
        END
