

opt=`cat config | awk 'BEGIN{FS="="}
/=/{ 
  opt[$1]=$2;
}
END{
  print opt["co_file"] " " opt["gr_file"] " " opt["iterations"];
}' `
./ARA $opt