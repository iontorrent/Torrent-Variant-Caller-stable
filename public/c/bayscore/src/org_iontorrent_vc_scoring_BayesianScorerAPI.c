 #include <jni.h>
 #include <stdio.h>
 #include "org_iontorrent_vc_scoring_BayesianScorerAPI.h"
 #include "posterior_flow.h"


 JNIEXPORT void JNICALL 
Java_org_iontorrent_vc_scoring_BayesianScorerAPI_rescorer_1init (JNIEnv *env, jobject o, jint x)
{
    rescorer_init(x);
}

JNIEXPORT jdouble JNICALL Java_org_iontorrent_vc_scoring_BayesianScorerAPI_calScore (JNIEnv *env, jobject o, jint pid)
{
     return calScore(pid);
}

JNIEXPORT void JNICALL Java_org_iontorrent_vc_scoring_BayesianScorerAPI_rescorer_1end (JNIEnv *env, jobject o)
{
    rescorer_end();
}

JNIEXPORT jint JNICALL Java_org_iontorrent_vc_scoring_BayesianScorerAPI_addRef (JNIEnv *env, jobject o, jint a, jbyteArray b)
{
     jboolean iscopy;
     jbyte *bb = env->GetByteArrayElements(b, &iscopy);
	return addRef(a, (char *) bb);
}
JNIEXPORT void JNICALL Java_org_iontorrent_vc_scoring_BayesianScorerAPI_addVariant__II_3B (JNIEnv *env, jobject o,jint pid, jint x, jbyteArray y)
{
     jboolean iscopy;
     jbyte *yy = env->GetByteArrayElements(y, &iscopy);
     addVariant(pid, x, (char *) yy);

}

JNIEXPORT void JNICALL Java_org_iontorrent_vc_scoring_BayesianScorerAPI_addVariant__II_3B_3B (JNIEnv *env, jobject o, jint pid, jint x, jbyteArray y, jbyteArray z)
{
     jboolean iscopy;
     jbyte *yy = env->GetByteArrayElements(y, &iscopy);
     jbyte *zz = env->GetByteArrayElements(z, &iscopy);
   	addVariant(pid, x, (char *) yy, (char *) zz); 
}

JNIEXPORT void JNICALL Java_org_iontorrent_vc_scoring_BayesianScorerAPI_addRead (JNIEnv *env, jobject o, jint pid, jint a, jint b, jint c, jbyteArray d, jintArray e, jbyteArray g, jint f)
{
     jboolean iscopy;
     jbyte *dd = env->GetByteArrayElements(d, &iscopy);
     jint *ee = env->GetIntArrayElements(e, &iscopy);
     jbyte *ff = env->GetByteArrayElements(g, &iscopy);

	addRead(pid, a, b, c, (char *) dd, (int *) ee, (char *) ff, f);
}

JNIEXPORT void JNICALL Java_org_iontorrent_vc_scoring_BayesianScorerAPI_finished (JNIEnv *env, jobject o, jint pid)
{
    finished(pid);
}
