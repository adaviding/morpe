﻿using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;

namespace Morpe.Validation
{
    public static class Chk
    {
        public static void Equal<T>(T observed, T expected, string message, params object[] args) where T : notnull
        {
            if (!expected.Equals(observed))
            {
                throw new ArgumentOutOfRangeException(string.Format(message, args));
            }
        }
        
        public static void Equal(object observed, object expected, string message, params object[] args)
        {
            if (observed == null && expected == null)
            {
                // do nothing
            }
            else if (observed == null
                     || expected == null
                     || !expected.Equals(observed))
            {
                throw new ArgumentOutOfRangeException(string.Format(message, args));
            }
        }
        
        public static void False(bool observed, string message, params object[] args)
        {
            if (observed)
            {
                throw new InvalidOperationException(string.Format(message, args));
            }
        }

        public static void Increasing<T>(IEnumerable<T> listing, string message, params object[] args)
            where T : System.IComparable<T>
        {
            T priorItem = default(T);
            bool firstTimeThrough = true;

            foreach (T item in listing)
            {
                if (firstTimeThrough)
                {
                    firstTimeThrough = false;
                }
                else
                {
                    if (priorItem.CompareTo(item) >= 0)
                    {
                        throw new InvalidOperationException(string.Format(message, args));
                    }
                }

                priorItem = item;
            }
        }

        public static void Less<T>(T lhs, T rhs, string message, params object[] args)
            where T : System.IComparable<T>
        {
            if ( 0 <= lhs.CompareTo( rhs ) )
            {
                throw new ArgumentOutOfRangeException(string.Format(message, args));
            }
        }
        
        public static void LessOrEqual<T>(T lhs, T rhs, string message, params object[] args)
            where T : System.IComparable<T>
        {
            if ( 0 < lhs.CompareTo( rhs ) )
            {
                throw new ArgumentOutOfRangeException(string.Format(message, args));
            }
        }
        
        public static void Finite(float observed, string message, params object[] args)
        {
            if (!float.IsFinite(observed))
            {
                throw new InvalidDataException(string.Format(message, args));
            }
        }
        
        public static void Finite(double observed, string message, params object[] args)
        {
            if (!double.IsFinite(observed))
            {
                throw new InvalidDataException(string.Format(message, args));
            }
        }

        public static void NotNull(object observed, string message, params object[] args)
        {
            if (observed == null)
            {
                throw new NullReferenceException(string.Format(message, args));
            }
        }

        public static void True(bool observed, string message, params object[] args)
        {
            if (!observed)
            {
                throw new InvalidOperationException(string.Format(message, args));
            }
        }
    }
}
