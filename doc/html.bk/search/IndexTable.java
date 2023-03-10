/*============================================================================

  Project: Simple JAVA Search Engine for Keyword Search
  
  JAVA Source file for the class IndexTable
  
  COPYRIGHT (C), 1998-2000, Thomas Baier, R Core Development Team

 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

  
  
  $Revision: 1.5 $

  $Date: 2000/02/10 17:03:55 $
  
  $Author: leisch $

============================================================================*/


import java.util.Vector;
import java.util.Enumeration;


public class IndexTable extends Vector
{

  public Vector search (String key, boolean searchDesc,
			boolean searchKeywords, boolean searchAliases)
  {
    Vector returnValue = new Vector ();
    Enumeration cursor = elements ();
    
    while (cursor.hasMoreElements ()) {
      IndexEntry entry = (IndexEntry) cursor.nextElement ();
      
      if (entry.matches (key, searchDesc,
			 searchKeywords, searchAliases)){
	returnValue.addElement (entry);
      }
    }
    if (!returnValue.isEmpty ()) {
      return returnValue;
    }
    
    return null;
  }
}


// Local Variables:
// mode: Java
// mode: font-lock
// End:
