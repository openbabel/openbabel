<xsl:stylesheet version = '1.0'
     xmlns:xsl='http://www.w3.org/1999/XSL/Transform'>
<xsl:output method="text"/>

<xsl:template match="/">
     <xsl:apply-templates select="//group"/>
</xsl:template>

<xsl:template match="group">
     <xsl:value-of select="@id"/>
     <xsl:text>&#10;</xsl:text>
     <xsl:value-of select="@Hall"/>
     <xsl:text>&#10;</xsl:text>
     <xsl:if test="@HMs">
          <xsl:value-of select="@HMs"/>
          <xsl:text>,</xsl:text>
     </xsl:if>
     <xsl:value-of select="@HM"/>
     <xsl:text>&#10;</xsl:text>
     <xsl:apply-templates select="transform"/>
     <xsl:text>&#10;</xsl:text>
</xsl:template>

<xsl:template match="transform">
     <xsl:value-of select="."/>
     <xsl:text>&#10;</xsl:text>
</xsl:template>

</xsl:stylesheet> 
