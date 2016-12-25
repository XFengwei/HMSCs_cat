source ~/.bash_profile
convert
cupid
cd ../../fitsDir/sedFitting/allSixMin/

rm *bg.fits

nLine=`awk 'END{print NR}' sourList.dat`
bands=''160' '250' '350' '500''
echo "this file has $nLine lines"
#
read -p "Please input the number of the starting line: " firstline

read -p "Please input the number of the ending line: " endline

for ((iline=$firstline; iline<=$endline; iline=iline+1))
	do
		sourName=`awk "{if (NR == $iline) print }" sourList.dat`
		
		rm *.sdf 

		for band in $bands
			do
				fits2ndf $sourName'_'$band'_smregrid36.fits' $sourName'_'$band'_raw.sdf' 
				findback $sourName'_'$band'_raw.sdf' $sourName'_'$band'_bg1.sdf' 14 RMS=! newalg=TRUE
				findback $sourName'_'$band'_bg1.sdf' $sourName'_'$band'_bg2.sdf' 14 RMS=! newalg=TRUE
				findback $sourName'_'$band'_bg2.sdf' $sourName'_'$band'_bg3.sdf' 14 RMS=! newalg=TRUE
				findback $sourName'_'$band'_bg3.sdf' $sourName'_'$band'_bg4.sdf' 14 RMS=! newalg=TRUE
				findback $sourName'_'$band'_bg4.sdf'  $sourName'_'$band'_bg.sdf' 14 RMS=! newalg=TRUE
				ndf2fits $sourName'_'$band'_bg.sdf' $sourName'_'$band'_bg.fits'
			done
		rm *.sdf

	done


