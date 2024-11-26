classdef ADRparameters
    properties
        D      % Diffusion coefficient
        v      % Velocity
        lam    % Lambda (decay rate or similar parameter)
        R      % Reaction rate or retardation factor
        gfun   % Function handle for a source term or boundary condition
        ffun   % Function handle for another function (e.g., sink or reaction)
        C0     % Initial concentration
    end

    methods
        % Constructor
        function obj = ADRparameters(D, v, lam, R, gfun, ffun, C0)
            if nargin > 0
                obj.D = D;
                obj.v = v;
                obj.lam = lam;
                obj.R = R;
                obj.gfun = gfun;
                obj.ffun = ffun;
                obj.C0 = C0;
            end
        end

        % Method to display parameters
        function displayParameters(obj)
            fprintf('D (Diffusion Coefficient): %f\n', obj.D);
            fprintf('v (Velocity): %f\n', obj.v);
            fprintf('lam (Lambda): %f\n', obj.lam);
            fprintf('R (Reaction Rate): %f\n', obj.R);
            fprintf('C0 (Initial Concentration): %f\n', obj.C0);
            fprintf('gfun: %s\n', func2str(obj.gfun));
            fprintf('ffun: %s\n', func2str(obj.ffun));
        end
    end
end
