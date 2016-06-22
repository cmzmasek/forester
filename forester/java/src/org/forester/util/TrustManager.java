package org.forester.util;

import java.io.IOException;
import java.net.URL;
import java.security.KeyManagementException;
import java.security.NoSuchAlgorithmException;
import java.security.cert.CertificateException;
import java.security.cert.X509Certificate;

import javax.net.ssl.HttpsURLConnection;
import javax.net.ssl.SSLContext;
import javax.net.ssl.X509TrustManager;

final class TrustManager implements X509TrustManager {
    
    @Override
    public X509Certificate[] getAcceptedIssuers() {
        return null;
    }

    @Override
    public void checkServerTrusted(X509Certificate[] paramArrayOfX509Certificate, String paramString)
        throws CertificateException {
    }

    @Override
    public void checkClientTrusted(X509Certificate[] paramArrayOfX509Certificate, String paramString)
        throws CertificateException {
    }
    
    final static HttpsURLConnection makeHttpsURLConnection( final URL url ) throws NoSuchAlgorithmException,
                                                                                        IOException,
                                                                                        KeyManagementException {
        
        final SSLContext ctx = SSLContext.getInstance("TLS");
        ctx.init(null, new TrustManager[] { new TrustManager() }, null);
        SSLContext.setDefault(ctx);
        
        final HttpsURLConnection con = (HttpsURLConnection) url.openConnection();
        return con;
        
    }
  
}
