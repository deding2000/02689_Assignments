classdef Flux
    properties
        a       % Parameter 'a'
        alpha   % Parameter 'alpha'
    end

    methods
        % Constructor
        function obj = Flux(a, alpha)
            if nargin > 0
                obj.a = a;
                obj.alpha = alpha;
            end
        end

        % Method to display parameters
        function displayParameters(obj)
            fprintf('a: %f\n', obj.a);
            fprintf('alpha: %f\n', obj.alpha);
        end
    end
end
