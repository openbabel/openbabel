/*
 *
 * FIPS-180-2 compliant SHA-256 implementation
 * Copyright (C) 2010 Paul Bakker.
 * Originally written (2003-2006) by Christophe Devine.
 *
 * This code is free software; you can redistribute it and/or
 * modify it, as a part of the International Chemical Identifier (InChI)
 * software, under the terms of the MIT license
 * https://opensource.org/license/mit
 *
 * Note that this licensing is only valid in combination with the InChI
 * software package, for all other cases contact Paul Bakker.
 */

/*
 *  The SHA-256 standard was published by NIST in 2002.
 *
 *  http://csrc.nist.gov/publications/fips/fips180-2/fips180-2.pdf
 */


#ifndef _SHA2_H_
#define _SHA2_H_

#ifdef __cplusplus
extern "C" {
#endif

/**
 * \brief          SHA-256 context structure
 */
    typedef struct
    {
        unsigned long total[2];     /*!< number of bytes processed  */
        unsigned long state[8];     /*!< intermediate digest state  */
        unsigned char buffer[64];   /*!< data block being processed */
    }
    sha2_context;

    /**
     * \brief          SHA-256 context setup
     *
     * \param ctx      SHA-256 context to be initialized
     */
    void sha2_starts( sha2_context *ctx );

    /**
     * \brief          SHA-256 process buffer
     *
     * \param ctx      SHA-256 context
     * \param input    buffer holding the  data
     * \param ilen     length of the input data
     */
    void sha2_update( sha2_context *ctx, unsigned char *input, int ilen );

    /**
     * \brief          SHA-256 final digest
     *
     * \param ctx      SHA-256 context
     * \param output   SHA-256 checksum result
     */
    void sha2_finish( sha2_context *ctx, unsigned char output[32] );

    /**
     * \brief          Output = SHA-256( input buffer )
     *
     * \param input    buffer holding the  data
     * \param ilen     length of the input data
     * \param output   SHA-256 checksum result
     */
    void sha2_csum( unsigned char *input, int ilen,
                    unsigned char output[32] );

    /**
     * \brief          Output = SHA-256( file contents )
     *
     * \param path     input file name
     * \param output   SHA-256 checksum result
     * \return         0 if successful, or 1 if fopen failed
     */
    int sha2_file( char *path, unsigned char output[32] );

    /**
     * \brief          Output = HMAC-SHA-256( input buffer, hmac key )
     *
     * \param key      HMAC secret key
     * \param keylen   length of the HMAC key
     * \param input    buffer holding the  data
     * \param ilen     length of the input data
     * \param output   HMAC-SHA-256 result
     */
    void sha2_hmac( unsigned char *key, int keylen,
                    unsigned char *input, int ilen,
                    unsigned char output[32] );

    /**
     * \brief          Checkup routine
     *
     * \return         0 if successful, or 1 if the test failed
     */
    int sha2_self_test( void );

#ifdef __cplusplus
}
#endif


#endif    /* _SHA2_H_ */

