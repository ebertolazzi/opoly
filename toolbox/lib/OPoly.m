classdef OPoly < handle

  properties (SetAccess = private, Hidden = true)
    objectHandle; % Handle to the underlying C++ class instance
  end

  methods

    function self = OPoly( varargin )
      % kind [, alpha, beta ]
      self.objectHandle = OPolyMexWrapper( 'new', varargin{:}  );
    end
    %
    % --------------------------------------------------------------------
    %
    function delete(self)
      % Destroy the C++ class instance
      OPolyMexWrapper( 'delete', self.objectHandle );
    end
    %
    % --------------------------------------------------------------------
    %
    function P = eval( self, n, x )
      P = OPolyMexWrapper( 'eval', self.objectHandle, n, x );
    end
    %
    % --------------------------------------------------------------------
    %
    function [P,Dp,s] = eval2( self, n, x )
      [P,Dp,s] = OPolyMexWrapper( 'eval2', self.objectHandle, n, x );
    end
    %
    % --------------------------------------------------------------------
    %
    function [PP,s] = sequence( self, n, x )
      [PP,s] = OPolyMexWrapper( 'sequence', self.objectHandle, n, x );
    end
    %
    % --------------------------------------------------------------------
    %
    function x = zeros( self, n, epsi )
      x = OPolyMexWrapper( 'zeros', self.objectHandle, n, epsi );
    end
    %
    % --------------------------------------------------------------------
    %
    function [node,w] = gauss( self, n, epsi )
      [node,w] = OPolyMexWrapper( 'gauss', self.objectHandle, n, epsi );
    end
  end
end
