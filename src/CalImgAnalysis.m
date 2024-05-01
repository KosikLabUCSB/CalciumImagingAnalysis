classdef CalImgAnalysis
% -------------------------------------------------------------------------
% CIA = CalImgAnalysis;
% -------------------------------------------------------------------------
% Ray - April 28, 2024
%
% Main class for Calcium Image Analysis for Kosik Lab.  This class contains 3 member variables 
% (also known as class properties): Preprocessing, 
% Classification, Visualization.  The data processing 
% functions are organized via these 3 member variables.
%
% To call a function, one would need to: first create a CalImgAnalysis object,
% then access the relevent member variable containing the function of 
% interest then finally make a call to the original function of interest.  
% For example, a call to the classification function findNeurons() 
% would look like this:
%
% CIA = CalImgAnalysis;
% CIA.Classification.findNeurons(X);

   properties
       Preprocessing;
       Classification;
       Visualization;
       Analysis;
   end
   methods
       function this = CalImgAnalysis()
           this.Preprocessing = Preprocessing;
           this.Classification = Classification;
           this.Visualization = Visualization;
           this.Analysis = Analysis;
       end
   end
end