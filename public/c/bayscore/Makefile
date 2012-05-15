all:	
	#rm -f src/org_iontorrent_vc_scoring_BayesianScorerAPI.h
	cd src; javac org/iontorrent/vc/scoring/BayesianScorerAPI.java
	cd src; javah -jni org.iontorrent.vc.scoring.BayesianScorerAPI
	#g++ -shared -o libBayesianScorerAPI.so -L/local/apps/jdk1.7.0_03/jre/lib/amd64/server -l jvm -fPIC posterior_flow.cpp util.cxx zutil.cxx fasta-io.cxx org_iontorrent_vc_scoring_BayesianScorerAPI.c -O2
	g++ -shared -o libBayesianScorerAPI.so -I/usr/lib/jvm/java-6-openjdk/include/ -L/usr/lib/jvm/java-6-openjdk/jre/lib/amd64/server/ -l jvm -fPIC src/posterior_flow.cpp src/util.cxx src/zutil.cxx src/fasta-io.cxx src/org_iontorrent_vc_scoring_BayesianScorerAPI.c -O2
	cd src; java -Djava.library.path=.. org.iontorrent.vc.scoring.BayesianScorerAPI
new:
	rm -f org_iontorrent_vc_scoring_BayesianScorerAPI.h
	javac org/iontorrent/vc/scoring/BayesianScorerAPI.java
	javah -jni org.iontorrent.vc.scoring.BayesianScorerAPI

