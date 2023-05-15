function vec = principalcomponent(domain,s,basis)

    switch(lower(basis))
        case 'laminar'
            x = domain.X;
            y = domain.Y;
            % Truncated (ky = -Q,0,Q) unstable solution to laminar state
            % See LSAofLaminarState.nb
            a1 = (-68/25-sqrt(323374)/25);
            a2 = (-2/5-sqrt(754)/5);
            a3 = (-4/27+sqrt(3166)/27);
            e1 = util.fftunstream(domain,a1.*cos(1.*x)+cos(1.*x-4.*y)-cos(1.*x+4.*y));
            e2 = util.fftunstream(domain,a2.*cos(2.*x)+cos(2.*x-4.*y)-cos(2.*x+4.*y));
            e3 = util.fftunstream(domain,a3.*cos(3.*x)+cos(3.*x-4.*y)-cos(3.*x+4.*y));
%             e1 = util.fftunstream(domain,cos(1.*x));
%             e2 = util.fftunstream(domain,cos(4.*y));
%             e3 = util.fftunstream(domain,cos(2.*x+4.*y));
            vec = [
                    util.fftinnerproduct(domain,e1,s);
                    util.fftinnerproduct(domain,e2,s);
                    util.fftinnerproduct(domain,e3,s);
                  ];
              
        case 'e0e2'
            x = domain.X;
            y = domain.Y;
            e1 = util.fftunstream(domain,cos(4.*y));
            e2 = util.fftunstream(domain,sin(2.*x));
            e3 = util.fftunstream(domain,sin(2.*x+4.*y));
            vec = [
                    util.fftinnerproduct(domain,e1,s);
                    util.fftinnerproduct(domain,e2,s);
                    util.fftinnerproduct(domain,e3,s);
                  ];
         case 'k1e1'
            x = domain.X;
            y = domain.Y;
%             e1 = util.load_eq(domain,'E1');
%             k1 = util.load_eq(domain,'K1_pace');
            v1 = util.fftunstream(domain,cos(4.*y));
            v2 = util.fftunstream(domain,sin(4.*y+1.*x));
            v3 = util.fftunstream(domain,sin(1.*y));
            vec = [
                    util.fftinnerproduct(domain,v1,s);
                    util.fftinnerproduct(domain,v2,s);
                    util.fftinnerproduct(domain,v3,s);
                  ];
        case 'e0e6'
            x = domain.X;
            y = domain.Y;
            e1 = util.fftunstream(domain,cos(4.*y));
            e2 = util.fftunstream(domain,sin(3.*x));
            e3 = util.fftunstream(domain,sin(3.*x+4.*y));
            vec = [
                    util.fftinnerproduct(domain,e1,s);
                    util.fftinnerproduct(domain,e2,s);
                    util.fftinnerproduct(domain,e3,s);
                  ];
        case 'e0e2e6'
            x = domain.X;
            y = domain.Y;
            e1 = util.fftunstream(domain,cos(4.*y));
            e2 = util.fftunstream(domain,sin(2.*x)+sin(3.*x));
            e3 = util.fftunstream(domain,sin(2.*x+4.*y)+sin(3.*x+4.*y));
            vec = [
                    util.fftinnerproduct(domain,e1,s);
                    util.fftinnerproduct(domain,e2,s);
                    util.fftinnerproduct(domain,e3,s);
                  ];
        case 'e10e04d'
            x = domain.X;
            y = domain.Y;
            e1 = util.fftunstream(domain,sin(0.*y)+sin(1.*x));
            e2 = util.fftunstream(domain,cos(0.*y)+cos(1.*x));
            vec = [
                    util.fftinnerproduct(domain,e1,s);
                    util.fftinnerproduct(domain,e2,s);
                    util.fftdissipation(domain,s);
                  ];
        case 'e0e2e6k1'
            x = domain.X;
            y = domain.Y;
            e1 = util.hat(domain,util.fftunstream(domain,cos(4.*y)));
            e2 = util.hat(domain,util.fftunstream(domain,sin(1.*x)+sin(2.*x)+sin(3.*x)));
            e3 = util.hat(domain,util.fftunstream(domain,sin(1.*x+4.*y)+sin(2.*x+4.*y)+sin(3.*x+4.*y)));
            vec = [
                    util.fftinnerproduct(domain,e1,s);
                    util.fftinnerproduct(domain,e2,s);
                    util.fftinnerproduct(domain,e3,s);
                  ];
        case 'poe1'
            obj = load('R40_PO01.mat','v1','v2','avgstate');
            vec = [
                    util.fftinnerproduct(domain,obj.v1,s);
                    util.fftinnerproduct(domain,obj.v2,s);
                    util.ffttwonorm(domain,s-obj.avgstate);
                  ];
        otherwise
            warning("basis not recognized or not supported: '%s'",basis)
    end

end