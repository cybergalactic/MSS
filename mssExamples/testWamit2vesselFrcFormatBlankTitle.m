function testWamit2vesselFrcFormatBlankTitle()
% testWamit2vesselFrcFormatBlankTitle() runs the FRC-layout check of
% testWamit2vesselFrcFormat on an alternative-2 file whose title line is
% blank and whose body mass is 30 (under 1000 in the file's units). The
% title line is still a line: frcFileAlternative must drop it and read the
% IMASS line, not the first row of the mass matrix.
%
% Usage: run from MATLAB with the MSS folders on the path.
%
% Author:    Enio Krizman
% Date:      2 November 2026

testWamit2vesselFrcFormat('', 30);

end
