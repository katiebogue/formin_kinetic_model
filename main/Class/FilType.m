classdef FilType
%FILTYPE Contains values for single, double, and dimer filament types
    %
    %   Construction:
    %       val=FILTYPE(single,double_a,dimer_a,double_b,dimer_b)
    %
    %       val=FILTYPE(single,double_a,dimer_a)
    %
    %   Multiple operators have been overloaded for this class (see methods).
    %
    %   See also FILAMENT.

    properties
        single % contains value for single filament type
        double % value for double filament type (filament class if 2)
        dimer  % value for dimer filament type (filament class if 2)
        intersectratio logical=false % if true, ratio will instead be overlapping values between double and dimer
    end
    properties (Dependent)
        ratio Filament % ratio of dimer/ double values (filament class; Dependent)
    end

    methods
        function obj = FilType(single,double_a,dimer_a,double_b,dimer_b)
            %FILTYPE Construct an instance of FilType
            %
            %   obj = FILTYPE(single,double,dimer)
            % 
            %   obj = FILTYPE(single,double_a,dimer_a,double_b,dimer_b)
            %
            %   assign single to single
            %   assign double to filament of double_a, double_b (or both
            %   double_a if no b provided)
            %   assign dimer to filament of dimer_a, dimer_b (or both
            %   dimer_a if no b provided)
            if nargin==5
                obj.single=single;
                obj.double=Filament(double_a,double_b);
                obj.dimer=Filament(dimer_a,dimer_b);
            elseif nargin==1
                obj.single=single;
                obj.double=single;
                obj.dimer=single;
            elseif nargin>0
                obj.single=single;
                obj.double=double_a;
                obj.dimer=dimer_a;
            end
        end

        function value = get.ratio(obj)
            % computes the ratio of the dimer/double objects, and does so
            % element by element
            value=getratio(obj);
            function out= getratio(input)
                if class(input.dimer)=="struct"
                    fieldsdimer=fieldnames(input.dimer);
                    fieldsdouble=fieldnames(input.double);
                    fields=intersect(fieldsdimer,fieldsdouble);
                    out=struct;
                    for i=1:length(fields)
                        temp1=FilType(0,input.double.(fields{i}),input.dimer.(fields{i}));
                        temp=getratio(temp1);                                                                                                                                                                                                                                
                        out.(fields{i})=temp;
                    end
                elseif class(input.dimer)=="Filament"
                    outa=getratio(FilType(0,input.double.a,input.dimer.a));
                    outb=getratio(FilType(0,input.double.b,input.dimer.b));
                    out=Filament(outa,outb);
                elseif class(input.dimer)=="string"
                    out="n/a";
                elseif length(input.double)>1 && length(input.double)==length(input.dimer)
                    if obj.intersectratio
                            out=intersect(input.dimer,input.double);
                    else
                        out=input.double;
                        for i=1:length(input.double)
                            fil=FilType(0,input.double(i),input.dimer(i));
                            out(i)=fil.ratio;
                        end
                    end
                else
                    if all(input.double==0) && all(input.dimer==0)
                        out=1;
                    elseif input.double==0
                        out=0;
                    elseif input.dimer==0
                        out=0;
                    else
                        if obj.intersectratio
                            out=intersect(input.dimer,input.double);
                        else
                            out=input.dimer./input.double;
                        end
                    end
                end
            end
        end

        function r = mrdivide(obj1,obj2)
            %MRDIVIDE divides obj1 by obj2 (matrix division)
            %
            %   filT=MRDIVIDE(fil1,fil2) generates FilType object by
            %   dividing (matrix division) fil1.single, fil1.double, and fil1.dimer by 
            %   fil2.single, fil2.double, and fil2.dimer,respectively
            %   Same as fil=fil1/fil2
            %
            %   filT=MRDIVIDE(fil,N) generates FilType object by dividing
            %   (matrix division) fil.single, fil.double, and fil.dimer by N
            %   Same as fil=fil/N
            %
            %   filT=MRDIVIDE(N,fil) generates FilType object by dividing 
            %   (matrix division) N by fil.single, fil.double, and fil.dimer
            %   Same as fil=N/fil
            %
            %   Inputs:
            %         obj1 : (double or FilType) the dividend; must be a 
            %                FilType if obj2 is a double
            %         obj2 : (double or FilType) the divisor; must be a 
            %                FilType if obj1 is a double
            %       
            %   Output is a FilType object.
            %
            %   See also FILTYPE, FILTYPE/USE_OPERATOR.
            operator=@mrdivide;
            r=use_operator(obj1,obj2,operator);
        end

        function r = plus(obj1,obj2)
            %PLUS adds obj1 and obj2
            %
            %   filT=PLUS(fil1,fil2) FilType filament object by
            %   adding fil1.single, fil1.double, and fil1.dimer with fil2.single, fil2.double, and fil2.dimer,
            %   respectively
            %   Same as fil=fil1+fil2
            %
            %   filT=PLUS(fil,N) generates FilType object by adding fil.single, fil.double, and fil.dimer with N
            %   Same as fil=fil+N
            %
            %   filT=PLUS(N,fil) generates FilType object by adding N with 
            %   fil.single, fil.double, and fil.dimer
            %   Same as fil=N+fil
            %
            %   Inputs:
            %         obj1 : (double or FilType) first object to be 
            %                added; must be a FilType if obj2 is a 
            %                double
            %         obj2 : (double or FilType) second object to be 
            %                added; must be a FilType if obj1 is a 
            %                double
            %       
            %   Output is a FilType object.
            %
            %   See also FILTYPE, FILTYPE/USE_OPERATOR.
            operator=@plus;
            r=use_operator(obj1,obj2,operator);
        end

        function r = minus(obj1,obj2)
            %MINUS subtracts obj2 from obj1
            %
            %   filT=MINUS(fil1,fil2) generates FilType object by 
            %   subtracting fil2.single, fil2.double, and fil2.dimer from fil1.single, fil1.double, and fil1.dimer,
            %   respectively
            %   Same as fil=fil1-fil2
            %
            %   filT=MINUS(fil,N) generates FilType object by subtracting 
            %   N from fil.single, fil.double, and fil.dimer
            %   Same as fil=fil-N
            %
            %   filT=MINUS(N,fil) generates FilType object by subtracting 
            %   fil.single, fil.double, and fil.dimer from N
            %   Same as fil=N-fil
            %
            %   Inputs:
            %         obj1 : (double or FilType) the base; must be 
            %                a FilType if obj2 is a double
            %         obj2 : (double or FilType) the exponent; must be
            %                a FilType if obj1 is a double
            %       
            %   Output is a FilType object.
            %
            %   See also FILTYPE, FILTYPE/USE_OPERATOR.
            operator=@minus;
            r=use_operator(obj1,obj2,operator);
        end

        function r = mtimes(obj1,obj2)
            %MTIMES applies matrix multiplication to obj1 and obj2
            %
            %   filT=MTIMES(fil1,fil2) generates FilType object by
            %   multiplying (matrix multiplication) fil1.single, fil1.double, and fil1.dimer with 
            %   fil2.single, fil2.double, and fil2.dimer,respectively
            %   Same as fil=fil1*fil2
            %
            %   filT=MTIMES(fil,N) generates FilType object by multiplying
            %   (matrix multiplication) fil.single, fil.double, and fil.dimer with N
            %   Same as fil=fil*N
            %
            %   filT=MTIMES(N,fil) generates FilType object by multiplying 
            %   (matrix multiplication) N with fil.single, fil.double, and fil.dimer
            %   Same as fil=N*fil
            %
            %   Inputs:
            %         obj1 : (double or FilType) first object to be 
            %                multiplied; must be a FilType if obj2 is a 
            %                double
            %         obj2 : (double or FilType) second object to be 
            %                multiplied; must be a FilType if obj1 is a 
            %                double
            %       
            %   Output is a FilType object.
            %
            %   See also FILTYPE, FILTYPE/USE_OPERATOR.
            operator=@mtimes;
            r=use_operator(obj1,obj2,operator);
        end
        

        function r = mpower(obj1,obj2)
            %MPOWER raises obj1 to the power of obj2
            %
            %   filT=MPOWER(fil1,fil2) generates FilType object by raising 
            %   fil1.single, fil1.double, and fil1.dimer to the powers of
            %   fil2.single, fil2.double, and fil2.dimer, respectively
            %   Same as fil=fil1^fil2
            %
            %   filT=MPOWER(fil,N) generates FilType object by raising 
            %   fil.single, fil.double, and fil.dimer to the power of N
            %   Same as fil=fil^N
            %
            %   filT=MPOWER(N,fil) generates FilType object by raising 
            %   N to the powers of fil.single, fil.double, and fil.dimer
            %   Same as fil=N^fil
            %
            %   Inputs:
            %         obj1 : (double or FilType) the base; must be 
            %                a FilType if obj2 is a double
            %         obj2 : (double or FilType) the exponent; must be
            %                a FilType if obj1 is a double
            %       
            %   Output is a FilType object.
            %
            %   See also FILTYPE, FILTYPE/USE_OPERATOR.
            operator=@mpower;
            r=use_operator(obj1,obj2,operator);
        end

        function r = exp(obj)
            %EXP applies the exponential (e^x) operator to a and b of a FILTYPE.
            %
            %   filT=EXP(obj) raises e to the power of obj and
            %   returns a new FILTYPE object
            %
            %   Inputs:
            %         obj : (FILTYPE) FILTYPE to apply exponential to
            %       
            %   Output is a FILTYPE object.
            %
            %   See also FILTYPE, FILTYPE/USE_OPERATOR.
            singler=exp(obj.single);
            doubler=exp(obj.double);
            dimerr=exp(obj.dimer);
            r=FilType(singler,doubler,dimerr);
        end


    end

    methods(Access = protected)
        function r = use_operator(obj1,obj2,operator)
            %USE_OPERATOR applies an operator to a FilType object and either another FilType object or a double
            %
            %   filT=USE_OPERATOR(obj1,obj2,operator) generates FilType
            %   object by applying operator to obj1 and obj2, element by
            %   element
            %
            %   Inputs:
            %         obj1     : (double or FilType) first object to apply
            %                    operation to; must be a FilType if obj2
            %                    is a double
            %         obj2     : (double or FilType) second object to apply
            %                    operation to; must be a FilType if obj1 
            %                    is a double
            %         operator : (function handle) operator to apply
            %       
            %   Output is a FilType object.
            % 
            %   See also FILTYPE/MRDIVIDE, FILTYPE/PLUS, FILTYPE/MINUS,
            %   FILTYPE/MTIMES, FILTYPE/MPOWER
            arguments
                obj1 
                obj2 
                operator function_handle
            end
            class1=class(obj1);
            class2=class(obj2);
            if class1=="double" && class2=="FilType"
                r=obj2;
                r.single=operator(obj1,obj2.single);
                if class(obj2.double)=="Filament"
                    r.double.a=operator(obj1,obj2.double.a);
                    r.double.b=operator(obj1,obj2.double.b);
                else
                    r.double=operator(obj1,obj2.double);
                end
                if class(obj2.dimer)=="Filament"
                    r.dimer.a=operator(obj1,obj2.dimer.a);
                    r.dimer.b=operator(obj1,obj2.dimer.b);
                else
                    r.dimer=operator(obj1,obj2.dimer);
                end
                return
            elseif class2=="FilType" && class1=="FilType"
                r=obj2;
                r.single=operator(obj1.single,obj2.single);
                if class(r.double)=="Filament" && class(obj1.double)=="Filament"
                    r.double.a=operator(obj1.double.a,obj2.double.a);
                    r.double.b=operator(obj1.double.b,obj2.double.b);
                else
                    r.double=operator(obj1.double,obj2.double);
                end
                if class(r.dimer)=="Filament" && class(obj1.dimer)=="Filament"
                    r.dimer.a=operator(obj1.dimer.a,obj2.dimer.a);
                    r.dimer.b=operator(obj1.dimer.b,obj2.dimer.b);
                else
                    r.dimer=operator(obj1.dimer,obj2.dimer);
                end
                return
            elseif class2=="double" && class1=="FilType"
                r=obj1;
                r.single=operator(obj1.single,obj2);
                if class(obj1.double)=="Filament"
                    r.double.a=operator(obj1.double.a,obj2);
                    r.double.b=operator(obj1.double.b,obj2);
                else
                    r.double=operator(obj1.double,obj2);
                end
                if class(obj1.dimer)=="Filament"
                    r.dimer.a=operator(obj1.dimer.a,obj2);
                    r.dimer.b=operator(obj1.dimer.b,obj2);
                else
                    r.dimer=operator(obj1.dimer,obj2);
                end
            else
                error("inputs must be of type FilType or double")
            end
        end
    end

    methods
        function obj=add_fils(obj)
            % ADD_FILS add a and b elements of double and dimer objects
            %
            %   dimer and double must both be of class Filament
            %   modifies the input object such that the double and dimer
            %   are no longer filaments but are the resulting sums.
            %
            % See also FILAMENT
            obj.double=obj.double.a+obj.double.b;
            obj.dimer=obj.dimer.a+obj.dimer.b;
        end

        function out=fil2array(obj)
            % FIL2ARRAY returns the FilType object as an array
            %
            %   Output is [single, double, dimer] or [single, double a,
            %   double b, dimer a, dimer b]
            %
            % See also FILAMENT
            if class(obj.dimer)=="Filament"
                dim=[obj.dimer.a,obj.dimer.b];
            else
                dim=obj.dimer;
            end
            if class(obj.double)=="Filament"
                dob=[obj.double.a,obj.double.b];
            else
                dob=obj.double;
            end
            out=[obj.single,dob,dim];
        end
    end
end