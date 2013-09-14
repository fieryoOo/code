@filelist=`ls XR_*.seed`;
$size=@filelist;
for ($i=0; $i<$size; $i++)
{
`rdseed -d -f $filelist[$i]\n`;
}
