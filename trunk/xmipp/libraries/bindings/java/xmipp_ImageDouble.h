/* DO NOT EDIT THIS FILE - it is machine generated */
#include <jni.h>
/* Header for class xmipp_ImageDouble */

#ifndef _Included_xmipp_ImageDouble
#define _Included_xmipp_ImageDouble
#ifdef __cplusplus
extern "C" {
#endif
/*
 * Class:     xmipp_ImageDouble
 * Method:    storeIds
 * Signature: ()V
 */
JNIEXPORT void JNICALL Java_xmipp_ImageDouble_storeIds
  (JNIEnv *, jclass);

/*
 * Class:     xmipp_ImageDouble
 * Method:    create
 * Signature: ()V
 */
JNIEXPORT void JNICALL Java_xmipp_ImageDouble_create
  (JNIEnv *, jobject);

/*
 * Class:     xmipp_ImageDouble
 * Method:    destroy
 * Signature: ()V
 */
JNIEXPORT void JNICALL Java_xmipp_ImageDouble_destroy
  (JNIEnv *, jobject);

/*
 * Class:     xmipp_ImageDouble
 * Method:    read
 * Signature: (Ljava/lang/String;)V
 */
JNIEXPORT void JNICALL Java_xmipp_ImageDouble_read
  (JNIEnv *, jobject, jstring);

/*
 * Class:     xmipp_ImageDouble
 * Method:    write
 * Signature: (Ljava/lang/String;)V
 */
JNIEXPORT void JNICALL Java_xmipp_ImageDouble_write
  (JNIEnv *, jobject, jstring);

/*
 * Class:     xmipp_ImageDouble
 * Method:    getData
 * Signature: ()[D
 */
JNIEXPORT jdoubleArray JNICALL Java_xmipp_ImageDouble_getData
  (JNIEnv *, jobject);

/*
 * Class:     xmipp_ImageDouble
 * Method:    getSize
 * Signature: ()[I
 */
JNIEXPORT jintArray JNICALL Java_xmipp_ImageDouble_getSize
  (JNIEnv *, jobject);

/*
 * Class:     xmipp_ImageDouble
 * Method:    printShape
 * Signature: ()V
 */
JNIEXPORT void JNICALL Java_xmipp_ImageDouble_printShape
  (JNIEnv *, jobject);

#ifdef __cplusplus
}
#endif
#endif
