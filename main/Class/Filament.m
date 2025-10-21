classdef Filament
%FILAMENT Contains values for a and b filament
    %
    %   Construction:
    %       fil=FILAMENT(a_val,b_val)
    %
    %   Multiple operators have been overloaded for this class (see methods).
    %
    %   See also FILTYPE.

    properties
        a  % value for filament a
        b  % value for filament b
    end
    properties(Dependent)
        avg % average of a and b values (Dependent)
    end

    methods
        function obj = Filament(a_val,b_val)
            %FILAMENT Construct an instance of Filament
            %
            %   fil=FILAMENT(a_val,b_val) generates filament object and 
            %       assigns a_val to a and b_val to b
            %
            %   fil=FILAMENT creates empty filament object
            %
            %   See also FILAMENT.
            if nargin>0
                obj.a = a_val;
                obj.b= b_val;
            end
        end

        function value = get.avg(obj)
            % computes average of obj.a and obj.b
            arr=[obj.a,obj.b];
            value=mean(arr);
        end

        function r = mrdivide(obj1,obj2)
            %MRDIVIDE divides obj1 by obj2 (matrix division)
            %
            %   fil=MRDIVIDE(fil1,fil2) generates filament object by
            %   dividing (matrix division) fil1.a and fil1.b by 
            %   fil2.a and fil2.b,respectively
            %   Same as fil=fil1/fil2
            %
            %   fil=MRDIVIDE(fil,N) generates filament object by dividing
            %   (matrix division) fil.a and fil.b by N
            %   Same as fil=fil/N
            %
            %   fil=MRDIVIDE(N,fil) generates filament object by dividing 
            %   (matrix division) N by fil.a and fil.b
            %   Same as fil=N/fil
            %
            %   Inputs:
            %         obj1 : (double or Filament) the dividend; must be a 
            %                Filament if obj2 is a double
            %         obj2 : (double or Filament) the divisor; must be a 
            %                Filament if obj1 is a double
            %       
            %   Output is a filament object.
            %
            %   See also FILAMENT, FILAMENT/USE_OPERATOR.
            operator=@mrdivide;
            r=use_operator(obj1,obj2,operator);
        end

        function r = plus(obj1,obj2)
            %PLUS adds obj1 and obj2
            %
            %   fil=PLUS(fil1,fil2) generates filament object by
            %   adding fil1.a and fil1.b with fil2.a and fil2.b,
            %   respectively
            %   Same as fil=fil1+fil2
            %
            %   fil=PLUS(fil,N) generates filament object by adding fil.a 
            %   and fil.b with N
            %   Same as fil=fil+N
            %
            %   fil=PLUS(N,fil) generates filament object by adding N with 
            %   fil.a and fil.b
            %   Same as fil=N+fil
            %
            %   Inputs:
            %         obj1 : (double or Filament) first object to be 
            %                added; must be a Filament if obj2 is a 
            %                double
            %         obj2 : (double or Filament) second object to be 
            %                added; must be a Filament if obj1 is a 
            %                double
            %       
            %   Output is a filament object.
            %
            %   See also FILAMENT, FILAMENT/USE_OPERATOR.
            operator=@plus;
            r=use_operator(obj1,obj2,operator);
        end

        function r = minus(obj1,obj2)
            %MINUS subtracts obj2 from obj1
            %
            %   fil=MINUS(fil1,fil2) generates filament object by 
            %   subtracting fil2.a and fil2.b from fil1.a and fil1.b,
            %   respectively
            %   Same as fil=fil1-fil2
            %
            %   fil=MINUS(fil,N) generates filament object by subtracting 
            %   N from fil.a and fil.b
            %   Same as fil=fil-N
            %
            %   fil=MINUS(N,fil) generates filament object by subtracting 
            %   fil.a and fil.b from N
            %   Same as fil=N-fil
            %
            %   Inputs:
            %         obj1 : (double or Filament) the base; must be 
            %                a Filament if obj2 is a double
            %         obj2 : (double or Filament) the exponent; must be
            %                a Filament if obj1 is a double
            %       
            %   Output is a filament object.
            %
            %   See also FILAMENT, FILAMENT/USE_OPERATOR.
            operator=@minus;
            r=use_operator(obj1,obj2,operator);
        end

        function r = mtimes(obj1,obj2)
            %MTIMES applies matrix multiplication to obj1 and obj2
            %
            %   fil=MTIMES(fil1,fil2) generates filament object by
            %   multiplying (matrix multiplication) fil1.a and fil1.b with 
            %   fil2.a and fil2.b,respectively
            %   Same as fil=fil1*fil2
            %
            %   fil=MTIMES(fil,N) generates filament object by multiplying
            %   (matrix multiplication) fil.a and fil.b with N
            %   Same as fil=fil*N
            %
            %   fil=MTIMES(N,fil) generates filament object by multiplying 
            %   (matrix multiplication) N with fil.a and fil.b
            %   Same as fil=N*fil
            %
            %   Inputs:
            %         obj1 : (double or Filament) first object to be 
            %                multiplied; must be a Filament if obj2 is a 
            %                double
            %         obj2 : (double or Filament) second object to be 
            %                multiplied; must be a Filament if obj1 is a 
            %                double
            %       
            %   Output is a filament object.
            %
            %   See also FILAMENT, FILAMENT/USE_OPERATOR.
            operator=@mtimes;
            r=use_operator(obj1,obj2,operator);
        end

        function r = mpower(obj1,obj2)
            %MPOWER raises obj1 to the power of obj2
            %
            %   fil=MPOWER(fil1,fil2) generates filament object by raising 
            %   fil1.a and fil1.b to the powers of fil2.a and fil2.b,
            %   respectively
            %   Same as fil=fil1^fil2
            %
            %   fil=MPOWER(fil,N) generates filament object by raising 
            %   fil.a and fil.b to the power of N
            %   Same as fil=fil^N
            %
            %   fil=MPOWER(N,fil) generates filament object by raising 
            %   N to the powers of fil.a and fil.b
            %   Same as fil=N^fil
            %
            %   Inputs:
            %         obj1 : (double or Filament) the base; must be 
            %                a Filament if obj2 is a double
            %         obj2 : (double or Filament) the exponent; must be
            %                a Filament if obj1 is a double
            %       
            %   Output is a filament object.
            %
            %   See also FILAMENT, FILAMENT/USE_OPERATOR.
            operator=@mpower;
            r=use_operator(obj1,obj2,operator);
        end

        function r = exp(obj)
            %EXP applies the exponential (e^x) operator to a and b of a filament.
            %
            %   fil=EXP(obj) raises e to the power of obj and
            %   returns a new filament object
            %
            %   Inputs:
            %         obj : (Filament) filament to apply exponential to
            %       
            %   Output is a filament object.
            %
            %   See also FILAMENT, FILAMENT/USE_OPERATOR.
            ar=exp(obj.a);
            br=exp(obj.b);
            r=Filament(ar,br);
        end
    end

    methods(Access = protected)
        function r = use_operator(obj1,obj2,operator)
            %USE_OPERATOR applies an operator to a Filament object and either another Filament object or a double
            %
            %   fil=USE_OPERATOR(obj1,obj2,operator) generates filament
            %   object by applying operator to obj1 and obj2
            %
            %   Inputs:
            %         obj1     : (double or Filament) first object to apply
            %                    operation to; must be a Filament if obj2
            %                    is a double
            %         obj2     : (double or Filament) second object to apply
            %                    operation to; must be a Filament if obj1 
            %                    is a double
            %         operator : (function handle) operator to apply
            %       
            %   Output is a filament object.
            %
            %   See also FILAMENT.
            arguments
                obj1 {mustBeA(obj1,["double","Filament"])}
                obj2 {mustBeA(obj2,["double","Filament"])}
                operator function_handle
            end
            if class(obj1)=="double"
                ar=operator(obj1,obj2.a);
                br=operator(obj1,obj2.b);
            elseif class(obj2)=="Filament"
                ar=operator(obj1.a,obj2.a);
                br=operator(obj1.b,obj2.b);
            elseif class(obj2)=="double"
                ar=operator(obj1.a,obj2);
                br=operator(obj1.b,obj2);
            end
            r=Filament(ar,br);
        end
    end
end